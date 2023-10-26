#!/usr/bin/env python
import logging
from intermol.exceptions import (UnimplementedFunctional, UnsupportedFunctional,
                                 UnimplementedSetting, UnsupportedSetting,
                                 GromacsError, InterMolError)

from collections import OrderedDict, defaultdict
logger = logging.getLogger('InterMolLog')
from intermol.system import System

class GromacsParser(object):
    class TopMoleculeType(object):
        """Inner class to store information about a molecule type."""
        def __init__(self):
            self.nrexcl = -1
            self.atoms = []
            self.bonds = []
            self.angles = []
            self.dihedrals = []
            self.rigidwaters = []
            self.exclusions = []
            self.pairs = []
            self.cmaps = []

    def __init__(self, top_file=None, gro_file=None,system=None):
        """Initializes the parser with all required metadata.

        Args:
            defines: Sets of default defines to use while parsing.
        """
        if not system:
            system = System()
        self.system = system

        self.top_filename = top_file
        self.gro_file = gro_file
        self.current_directive = None
        self.if_stack = list()
        self.else_stack = list()
        self.molecule_types = OrderedDict()
        self.molecules = list()
        self.current_molecule_type = None
        self.current_molecule = None
        self.bondtypes = dict()
        self.angletypes = dict()
        self.dihedraltypes = dict()
        self.implicittypes = dict()
        self.pairtypes = dict()
        self.cmaptypes = dict()
        self.nonbondedtypes = dict()

    def too_few_fields(self, line):
        raise GromacsError('Too few fields in [ {0} ] line: {1}'.format(
            self.current_directive, line))

    def invalid_line(self, line):
        raise GromacsError('Invalid format in [ {0} ] line: {1}'.format(
            self.current_directive, line))

    # =========== Pre-processing and forcetype creation =========== #
    def process_file(self, top_filename):
        append = ''
        with open(top_filename) as top_file:
            for line in top_file:
                if line.strip().endswith('\\'):
                    append = '{0} {1}'.format(append, line[:line.rfind('\\')])
                else:
                    self.process_line(top_filename, '{0} {1}'.format(append, line))
                    append = ''

    def process_line(self, top_filename, line):
        """Process one line from a file."""
        if ';' in line:
            line = line[:line.index(';')]
        stripped = line.strip()
        ignore = not all(self.if_stack)
        if stripped.startswith('*') or len(stripped) == 0:
            # A comment or empty line.
            return

        elif stripped.startswith('[') and not ignore:
            # The start of a category.
            if not stripped.endswith(']'):
                raise GromacsError('Illegal line in .top file: '+line)
            self.current_directive = stripped[1:-1].strip()
            logger.debug("Parsing {0}...".format(self.current_directive))

        elif not ignore:
            # A line of data for the current category
            if self.current_directive is None:
                raise GromacsError('Unexpected line in .top file: "{0}"'.format(line))
            if self.current_directive == 'defaults':
                self.process_defaults(line)
            elif self.current_directive == 'moleculetype':
                self.process_moleculetype(line)
            elif self.current_directive == 'molecules':
                self.process_molecule(line)
            elif self.current_directive == 'atoms':
                self.process_atom(line)
            elif self.current_directive == 'bonds':
                self.process_bond(line)
            elif self.current_directive == 'angles':
                self.process_angle(line)
            elif self.current_directive == 'dihedrals':
                self.process_dihedral(line)
            elif self.current_directive == 'settles':
                self.process_settle(line)
            elif self.current_directive == 'exclusions':
                self.process_exclusion(line)
            elif self.current_directive == 'pairs':
                self.process_pair(line)
            elif self.current_directive == 'cmap':
                self.process_cmap(line)
            elif self.current_directive == 'atomtypes':
                self.process_atomtype(line)
            elif self.current_directive == 'bondtypes':
                self.process_bondtype(line)
            elif self.current_directive == 'angletypes':
                self.process_angletype(line)
            elif self.current_directive == 'dihedraltypes':
                self.process_dihedraltype(line)
            elif self.current_directive == 'implicit_genborn_params':
                self.process_implicittype(line)
            elif self.current_directive == 'pairtypes':# and not self.system.genpairs:
                self.process_pairtype(line)
            elif self.current_directive == 'cmaptypes':
                self.process_cmaptype(line)
            elif self.current_directive == 'nonbond_params':
                self.process_nonbond_params(line)
            elif self.current_directive.startswith('virtual_sites'):
                vsite_type = self.current_directive[-1]
                self.process_virtual_sites(line, vsite_type)

    def process_defaults(self, line):
        """Process the [ defaults ] line."""
        fields = line.split()
        if len(fields) < 3:
            self.too_few_fields(line)
        self.system.nonbonded_function = int(fields[0])
        self.system.combination_rule = int(fields[1])
        self.system.genpairs = fields[2]
        #self.system.lj_correction = float(fields[3])
        #self.system.coulomb_correction = float(fields[4])

    def process_moleculetype(self, line):
        """Process a line in the [ moleculetypes ] category."""
        fields = line.split()
        if len(fields) < 1:
            self.too_few_fields(line)
        mol_type = self.TopMoleculeType()
        mol_type.nrexcl = int(fields[1])
        self.molecule_types[fields[0]] = mol_type
        self.current_molecule_type = mol_type

    def process_molecule(self, line):
        """Process a line in the [ molecules ] category."""
        fields = line.split()
        if len(fields) < 2:
            self.too_few_fields(line)
        self.molecules.append((fields[0], int(fields[1])))

    def process_atom(self, line):
        """Process a line in the [ atoms ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)
        if len(fields) not in [6, 7, 8, 11]:
            self.invalid_line(line)
        self.current_molecule_type.atoms.append(fields)

    def process_bond(self, line):
        """Process a line in the [ bonds ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 3:
            self.too_few_fields(line)
        self.current_molecule_type.bonds.append(fields)

    def process_angle(self, line):
        """Process a line in the [ angles ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 4:
            self.too_few_fields(line)
        self.current_molecule_type.angles.append(fields)

    def process_dihedral(self, line):
        """Process a line in the [ dihedrals ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)
        self.current_molecule_type.dihedrals.append(fields)

    def process_settle(self, line):
        """Process a line in the [ settles ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 4:
            self.too_few_fields(line)
        self.current_molecule_type.rigidwaters.append(fields)

    def process_exclusion(self, line):
        """Process a line in the [ exclusions ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 2:
            self.too_few_fields(line)
        self.current_molecule_type.exclusions.append(fields)

    def process_pair(self, line):
        """Process a line in the [ pairs ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype()
        fields = line.split()
        if len(fields) < 3:
            self.too_few_fields(line)
        self.current_molecule_type.pairs.append(fields)

    def process_cmap(self, line):
        """Process a line in the [ cmaps ] category."""
        if self.current_molecule_type is None:
            self.directive_before_moleculetype('cmap')
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        self.current_molecule_type.cmaps.append(fields)

    def process_atomtype(self, line):
        """Process a line in the [ atomtypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        if len(fields[3]) == 1:
            # Bonded type and atomic number are both missing.
            fields.insert(1, None)
            fields.insert(1, None)
        elif len(fields[4]) == 1 and len(fields[5]) >= 1:
            if fields[1][0].isalpha():
                # Atomic number is missing.
                fields.insert(2, None)
            else:
                # Bonded type is missing.
                fields.insert(1, None)

        atomtype = fields[0]
        if fields[1] == None:
            bondingtype = atomtype
        else:
            bondingtype = fields[1]
        if fields[2]:
            atomic_number = int(fields[2])
        else:
            atomic_number = -1
        mass = float(fields[3])  #* units.amu
        charge = float(fields[4]) # * units.elementary_charge
        ptype = fields[5]
        ## Add correct units to the LJ parameters.
        #if self.system.combination_rule == "Multiply-C6C12":
            #lj_param1 = (float(fields[6]) *
                         #units.kilojoules_per_mole * units.nanometers**(6))
            #lj_param2 = (float(fields[7]) *
                         #units.kilojoules_per_mole * units.nanometers**(12))
            #AtomtypeClass = AtomCType
        #elif self.system.combination_rule in ['Multiply-Sigeps', 'Lorentz-Berthelot']:
            #lj_param1 = float(fields[6]) * units.nanometers           # sigma
            #lj_param2 = float(fields[7]) * units.kilojoules_per_mole  # epsilon
            #AtomtypeClass = AtomSigepsType
        #else:
            #raise InterMolError("Unknown combination rule: {0}".format(self.system.combination_rule))
        #new_atom_type = AtomtypeClass(atomtype, bondingtype, atomic_number,
                                      #mass, charge, ptype, lj_param1, lj_param2)
        #self.system.add_atomtype(new_atom_type)

    def process_bondtype(self, line):
        """Process a line in the [ bondtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)

        btypes = fields[:2]
        bond_type = self.process_forcetype(btypes, 'bond', line, 2,
                self.gromacs_bond_types, self.canonical_bond)
        self.bondtypes[tuple(fields[:2])] = bond_type

    def process_angletype(self, line):
        """Process a line in the [ angletypes ] category."""
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        btypes = fields[:3]
        angle_type = self.process_forcetype(btypes, 'angle', line, 3,
                self.gromacs_angle_types, self.canonical_angle)
        self.angletypes[tuple(fields[:3])] = angle_type

    def process_dihedraltype(self, line):
        """Process a line in the [ dihedraltypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)

        # Some gromacs parameters don't include sufficient numbers of types.
        # Add some zeros (bit of a kludge).
        line += ' 0.0 0.0 0.0'
        fields = line.split()

        # Check whether they are using 2 or 4 atom types
        if fields[2].isdigit():
            btypes = ['X', fields[0], fields[1], 'X']
            n_atoms_specified = 2
        elif fields[4].isdigit() and not fields[3].isdigit(): # assumes gromacs types are not all digits.
            btypes = fields[:4]
            n_atoms_specified = 4
        else:
            # TODO: Come up with remaining cases (are there any?) and a proper
            #       failure case.
            logger.warning('Should never have gotten here.')
        dihedral_type = self.process_forcetype(
            btypes, 'dihedral', line, n_atoms_specified,
            self.gromacs_dihedral_types, self.canonical_dihedral)

        # Still need a bit more information
        numeric_dihedraltype = fields[n_atoms_specified]
        dihedral_type.improper = numeric_dihedraltype in ['2', '4']

        key = tuple([btypes[0], btypes[1], btypes[2], btypes[3],
                     dihedral_type.improper])

        if key in self.dihedraltypes:
            # There are multiple dihedrals defined for these atom types.
            self.dihedraltypes[key].add(dihedral_type)
        else:
            self.dihedraltypes[key] = {dihedral_type}

    def process_forcetype(self, bondingtypes, forcename, line, n_atoms,
                          gromacs_force_types, canonical_force):
        """ """
        fields = line.split()

        numeric_forcetype = fields[n_atoms]
        gromacs_force_type = gromacs_force_types[numeric_forcetype]
        kwds = self.create_kwds_from_entries(fields, gromacs_force_type, offset=n_atoms+1)
        CanonicalForceType, kwds = canonical_force(
            kwds, gromacs_force_type, direction='into')

        force_type = CanonicalForceType(*bondingtypes, **kwds)

        if not force_type:
            logger.warning("{0} is not a supported {1} type".format(fields[2], forcename))
            return
        else:
            return force_type

    def process_implicittype(self, line):
        """Process a line in the [ implicit_genborn_params ] category."""
        fields = line.split()
        if len(fields) < 6:
            self.too_few_fields(line)
        self.implicittypes[fields[0]] = fields

    def process_pairtype(self, line):
        """Process a line in the [ pairtypes ] category."""
        fields = line.split()
        if len(fields) < 5:
            self.too_few_fields(line)

        pair_type = None
        PairFunc = None
        combination_rule = self.system.combination_rule
        kwds = dict()
        numeric_pairtype = fields[2]
        if numeric_pairtype == '1':
            # LJ/Coul. 1-4 (Type 1)
            if len(fields) == 5:
                if combination_rule == "Multiply-C6C12":
                    PairFunc = LjCPairType
                elif combination_rule in ['Multiply-Sigeps', 'Lorentz-Berthelot']:
                    PairFunc = LjSigepsPairType
            offset = 3
        elif numeric_pairtype == '2':
            if combination_rule == "Multiply-C6C12":
                PairFunc = LjqCPairType
            elif combination_rule in ['Multiply-Sigeps', 'Lorentz-Berthelot']:
                PairFunc = LjqSigepsPairType
            offset = 4
        else:
            logger.warning("Could not find pair type for line: {0}".format(line))

        if PairFunc:
            pairvars = [fields[0], fields[1]]
            kwds = self.create_kwds_from_entries(fields, PairFunc, offset=offset)
            # kludge because of placement of scaleQQ...
            if numeric_pairtype == '2':
                # try to get this out ...
                kwds['scaleQQ'] = float(fields[3]) * units.dimensionless
            pair_type = PairFunc(*pairvars, **kwds)

        self.pairtypes[tuple(fields[:2])] = pair_type

    def process_cmaptype(self, line):
        """Process a line in the [ cmaptypes ] category."""
        fields = line.split()
        if len(fields) < 8 or len(fields) < 8+int(fields[6])*int(fields[7]):
            self.too_few_fields(line)
        self.cmaptypes[tuple(fields[:5])] = fields

    def process_nonbond_params(self, line):
        """Process a line in the [ nonbond_param ] category."""
        fields = line.split()
        NonbondedFunc = None
        combination_rule = self.system.combination_rule

        if fields[2] == '1':
            if combination_rule == 'Multiply-C6C12':
                NonbondedFunc = LjCNonbondedType
            elif combination_rule in ['Lorentz-Berthelot', 'Multiply-Sigeps']:
                NonbondedFunc = LjSigepsNonbondedType
        elif fields[2] == '2':
            if combination_rule == 'Buckingham':
                NonbondedFunc = BuckinghamNonbondedType
        else:
            logger.warning("Could not find nonbonded type for line: {0}".format(line))

        nonbonded_vars = [fields[0], fields[1]]
        kwds = self.create_kwds_from_entries(fields, NonbondedFunc, offset=3)
        nonbonded_type = NonbondedFunc(*nonbonded_vars, **kwds)
        # TODO: figure out what to do with the gromacs numeric type
        nonbonded_type.type = int(fields[2])
        self.system.nonbonded_types[tuple(nonbonded_vars)] = nonbonded_type

    def process_virtual_sites(self, line, v_site_type):
        """Process a line in a [ virtual_sites? ] category."""
        if v_site_type == 'n':
            raise UnimplementedSetting('Parsing of [ virtual_sitesn ] directives'
                                       ' is not yet implemented')
        fields = line.split()
        self.current_molecule_type.virtuals[v_site_type].append(fields)
