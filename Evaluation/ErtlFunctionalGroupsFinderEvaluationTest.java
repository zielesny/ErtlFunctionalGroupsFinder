/**
 * Evaluation tests for
 * ErtlFunctionalGroupsFinder for CDK
 * Copyright (C) 2019 Jonas Schaub
 * 
 * Source code is available at <https://github.com/zielesny/ErtlFunctionalGroupsFinder>
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.openscience.cdk.tools;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.hash.AtomEncoder;
import org.openscience.cdk.hash.BasicAtomEncoder;
import org.openscience.cdk.hash.HashGeneratorMaker;
import org.openscience.cdk.hash.MoleculeHashGenerator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ErtlFunctionalGroupsFinder.Mode;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;


/**
 * This test class can be used to read an SD file containing chemical structures, to extract their functional groups using
 * the ErtlFunctionalGroupsFinder with different settings (i.e. electron donation model and cycle finder algorithm), and write 
 * the functional groups with their associated frequency under the given settings in this SD file to a CSV file.
 * <p>
 * To run correctly the constant SD_FILE_PATH must be set to where to find the specific file on the local system.
 * <p>
 * All written files will be placed in a new folder in the same directory as the read SD file.
 * <p>
 * Note for addition of new tests: Only one SD file should be analyzed per test method (since some mechanisms work under 
 * that assumption).
 * 
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class ErtlFunctionalGroupsFinderEvaluationTest extends CDKTestCase {
    
    //<editor-fold defaultstate="collapsed" desc="Private static final constants">
    
    //<editor-fold defaultstate="collapsed" desc="Paths, directories and files">
    /**
     * Path to SD file that should be analyzed
     */
    private static final String SD_FILE_PATH = "...\\ChEBI_lite_3star_subset.sdf";
    
    /**
     * Directory for output files; Will be created as sub-folder in the working directory (the directory of the read SD file)
     */
    private static final String OUTPUT_FOLDER_FROM_WORKING_DIRECTORY = "ErtlFunctionalGroupsFinderEvaluationTest_Output";
    
    /**
     * Format of the time stamp addition to all written output files
     */
    private static final String DATE_TIME_FORMAT_PATTERN = "uuuu_MM_dd_HH_mm";
    
    /**
     * Separator for file name segments (test identifier, file name, time stamp)
     */
    private static final String FILE_NAME_ADDITION_SEPERATOR = "_";
    
    //<editor-fold defaultstate="collapsed" desc="Test identifiers">
    /**
     * Identifier string for SD file electron donation test that counts multiples
     */
    private static final String ELECTRON_DONATION_TEST_IDENTIFIER = "ElectronDonationTest";
    
    /**
     * Identifier string for SD file electron donation test that does not count multiples per molecule
     */
    private static final String ELECTRON_DONATION_NO_MULTIPLES_TEST_IDENTIFIER = "ElectrondDonationNoMultiplesTest";
    
    /**
     * Identifier string for test of different CycleFinder settings
     */
    private static final String CYCLE_FINDER_TEST_IDENTIFIER = "CycleFinderTest";
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Exceptions log file">
    /**
     * Name of file for logging occurred exceptions with causing molecules
     */
    private static final String EXCEPTIONS_LOG_FILE_NAME = "Exceptions_Log";
    
    /**
     * File type of exceptions log file
     */
    private static final String EXCEPTIONS_LOG_FILE_TYPE = ".txt";
    
    /**
     * First lines in the exceptions log file
     */
    private static final String EXCEPTIONS_LOG_FILE_HEADER = "Following molecules led to the specified exceptions:"
            + System.getProperty("line.separator")
            + "(Note: If too many exceptions are thrown too fast the JVM stops filling in the complete stack trace. "
            + "You need to be looking at an earlier stack trace to see the details.)"
            + System.getProperty("line.separator");
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Filtered molecules log file">
    /**
     * Name of file for logging filtered molecules
     */
    private static final String FILTERED_MOLECULES_FILE_NAME = "Filtered_Molecules";
    
    /**
     * File type of filtered molecules log file
     */
    private static final String FILTERED_MOLECULES_FILE_TYPE = ".txt";
    
    /**
     * First line in the filtered molecules log file
     */
    private static final String FILTERED_MOLECULES_FILE_HEADER = "Following molecules were filtered:\n";
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Output file">
    /**
     * Name of file for writing resulting functional groups and frequencies to (output file)
     */
    private static final String OUTPUT_FILE_NAME = "Functional_Groups_Frequencies";
    
    /**
     * File type of output file
     */
    private static final String OUTPUT_FILE_TYPE = ".csv";
    
    /**
     * Key for the output file's header under which to store the unique SMILES code of a functional group
     */
    private static final String SMILES_CODE_KEY = "smiles";
    
    /**
     * Key for the output file's header under which to store the pseudo SMILES code of a functional group
     */
    private static final String PSEUDO_SMILES_CODE_KEY = "pseudoSmiles";
    
    /**
     * Key for the output file's header under which to store the hash code of a functional group
     */
    private static final String HASH_CODE_KEY = "hashCode";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * cdk electron donation model (and internally for the master HashMap's inner maps)
     */
    private static final String CDK_MODEL_SETTINGS_KEY = "cdk";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * daylight electron donation model (and internally for the master HashMap's inner maps)
     */
    private static final String DAYLIGHT_MODEL_SETTINGS_KEY = "daylight";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * pi bonds electron donation model (and internally for the master HashMap's inner maps)
     */
    private static final String PIBONDS_MODEL_SETTINGS_KEY = "piBonds";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * cdk allowing exocyclic electron donation model (and internally for the master HashMap's inner maps)
     */
    private static final String CDK_EXOCYCLIC_MODEL_SETTINGS_KEY = "cdkExocyclic";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * cdk legacy aromaticity model (and internally for the master HashMap's inner maps)
     */
    private static final String CDK_LEGACY_MODEL_SETTINGS_KEY = "cdkLegacy";
    
    /**
     * This string will be added to an original settings key when using the ErtlFunctionalGroupsFinder with activated 
     * generalization
     */
    private static final String GENERALIZATION_SETTINGS_KEY_ADDITION = "Generalized";
    
    /**
     * Key for the master HashMap's inner maps under which to store the ChEBI or ChEMBL id/name or CDK title of an exemplary 
     * molecule that contains this functional group
     */
    private static final String MOLECULE_OF_ORIGIN_KEY = "origin";
    
    /**
     * Separator for the output file's values
     */
    private static final String OUTPUT_FILE_SEPERATOR = ",";
    
    /**
     * Placeholder String for every functional group's SMILES code whose real SMILES representation could not be 
     * generated
     */
    private static final String SMILES_CODE_PLACEHOLDER = "[SMILES code could not be created]";
    
    /**
     * Placeholder String for 'parent' molecules' IDs, assigned in case they could not be read
     */
    private static final String MOLECULE_OF_ORIGIN_ID_PLACEHOLDER = "[no id for molecule of origin]";
    //</editor-fold>
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Filtering and Preprocessing">
    /**
     * Key for setting an IAtomContainer's property that specifies if the molecule consisted of one or more unconnected
     * structures and the biggest of these structures was selected in preprocessing
     */
    private static final String BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY = "biggest_fragment_was_selected";
    
    /**
     * Key for setting an IAtomContainer's property that specifies if the molecule contained charged atoms and these
     * charges were neutralized in preprocessing
     */
    private static final String CHARGES_NEUTRALIZED_PROPERTY_KEY = "charges_were_neutralized";
    
    /**
     * Key for setting an IAtomContainers property that specifies if the molecule must be filtered
     */
    private static final String MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY = "molecule_must_be_filtered";
    
    /**
     * Key for setting an IAtomContainers property that specifies why the molecule was or has to be filtered
     */
    private static final String CAUSE_FOR_FILTERING_PROPERTY_KEY = "filtering_cause";
    
    /**
     * Message specifying that the atom or bond count of a molecule is zero
     */
    private static final String ATOM_OR_BOND_COUNT_ZERO = "Atom or bond count 0";
    
    /**
     * Message specifying that a molecule contains not allowed atomic numbers
     */
    private static final String FORBIDDEN_ATOMIC_NUMBER = "Contains one or more metal, metalloid or \"R\" atoms";
    
    /**
     * All allowed atomic numbers to pass to the ErtlFunctionalGroupsFinder; 
     * String will be split and resulting integers passed to a set
     */
    private static final String NON_METALLIC_ATOMIC_NUMBERS = "1,2,6,7,8,9,10,15,16,17,18,34,35,36,53,54,86";
    
    /**
     * This string will be added to an original settings key (only for exception logging) when a molecule consists of two or 
     * more unconnected structures and the biggest one is chosen for analysis in the preprocessing
     */
    private static final String FRAGMENT_SELECTED_SETTINGS_KEY_ADDITION = "(biggest fragment selected)";
    
    /**
     * This string will be added to an original settings key (only for exception logging) when a molecule contains charged atoms
     * and these charges are neutralized in the preprocessing
     */
    private static final String NEUTRALIZED_SETTINGS_KEY_ADDITION = "(neutralized)";
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Pseudo SMILES code">
    /**
     * Pseudo SMILES representation of an aromatic C atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_CARBON = "C*";
    
    /**
     * Pseudo SMILES representation of an aromatic N atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_NITROGEN = "N*";
    
    /**
     * Pseudo SMILES representation of an aromatic S atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_SULPHUR = "S*";
    
    /**
     * Pseudo SMILES representation of an aromatic O atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_OXYGEN = "O*";
    
    /**
     * Pseudo SMILES representation of an aromatic Se atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_SELENIUM = "Se*";
    
    /**
     * Pseudo SMILES representation of an aromatic P atom
     */
    private static final String PSEUDO_SMILES_AROMATIC_PHOSPHOR = "P*";
    
    /**
     * Pseudo SMILES representation of an undefined pseudo atom
     */
    private static final String PSEUDO_SMILES_R_ATOM = "R";
    //</editor-fold>
    
    /**
     * Type of the generated SMILES codes
     */
    private static final int SMILES_GENERATOR_OUTPUT_MODE = SmiFlavor.Unique;
    
    /**
     * Initial capacity of the master HashMap (where all data is written to); May be adjusted when analyzing larger 
     * SD files
     */
    private static final int MASTER_HASHMAP_INITIAL_CAPACITY = 1000;
    
    /**
     * Initial capacity of the master HashMap's inner maps that store the frequencies for different settings for a 
     * specific functional group
     */
    private static final int INNER_HASHMAPS_INITIAL_CAPACITY = 20;
    
    /**
     * Load factor of the master HashMap
     */
    private static final float MASTER_HASHMAP_LOAD_FACTOR = 0.9f;
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Directory for all produced files; It will be the directory where th SD file that is analyzed was loaded from
     */
    private String outputDirectory;
    
    /**
     * Counts all occurring exceptions in one test
     */
    private int exceptionsCounter;
    
    /**
     * True if the filtered molecules were logged in the filtered molecules log file; This is only necessary in the 
     * first iteration since the applied filters are the same in every iteration (assuming that in a single test
     * only one SD file is analyzed)
     */
    private boolean areFilteredMoleculesLogged;
    
    /**
     * True if all operations in initialize() were successful and the test is able to run
     */
    private boolean isTestAbleToRun;
    
    /**
     * True if the currently running test requires files to be loaded and files to log exceptions or filtered molecules 
     * or the produced results and therefore the PrintWriter objects have been initialized
     */
    private boolean areFileOperationsActivated;
    
    /**
     * PrintWriter for logging exceptions
     */
    private PrintWriter exceptionsPrintWriter;
    
    /**
     * PrintWriter for logging filtered molecules
     */
    private PrintWriter filteredMoleculesPrintWriter;
    
    /**
     * PrintWriter for writing the resulting functional groups and frequencies
     */
    private PrintWriter dataOutputPrintWriter;
    
    /**
     * SmilesGenerator for generating SMILES and pseudo SMILES representations
     */
    private SmilesGenerator smilesGenerator;
    
    /**
     * MoleculeHashGenerator for easy assessment whether a functional group was already entered into the master HashMap
     */
    private MoleculeHashGenerator molHashGenerator;
    
    /**
     * Instance of the ErtlFunctionalGroupsFinder with generalization turned off
     */
    private ErtlFunctionalGroupsFinder ertlFGFinderGenOff;
    
    /**
     * Instance of the ErtlFunctionalGroupsFinder with generalization turned on
     */
    private ErtlFunctionalGroupsFinder ertlFGFinderGenOn;
    
    /**
     * Master HashMap for storing results; Its keys are the hash codes produced by the MoleculeHashGenerator for the 
     * functional groups and its values are inner HashMaps that hold the (pseudo) SMILES representation of a functional 
     * group and its frequencies for different settings as String-Object pairs, plus an exemplary molecule of origin
     */
    private HashMap<Long, HashMap> masterHashMap;
    
    /**
     * A map that gives a certain element symbol for a placeholder atom marking a specific aromatic atom in pseudo SMILES
     * creation
     */
    private HashMap<String, String> pseudoSmilesAromaticElementToPlaceholderElementMap;
    
    /**
     * A map that gives the pseudo SMILES representation for a specific placeholder element from 
     * pseudoSmilesAromaticElementToPlaceholderElementMap
     */
    private HashMap<String, String> pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap;
    
    /**
     * A list for storing all used settings keys in a test
     */
    private List<String> settingsKeysList;
    
    /**
     * All allowed atomic numbers to pass to the ErtlFunctionalGroupsFinder as a set of integers (will be parsed from 
     * NON_METALLIC_ATOMIC_NUMBERS)
     */
    private Set<Integer> nonMetallicAtomicNumbersSet;
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * <p>
     * Note: it does not initialize any class variables (except 5) because that would be unnecessary when it is called by a 
     * test method inherited from CDKTestCase; these initializations are done by initialize().
     */
    public ErtlFunctionalGroupsFinderEvaluationTest() {
        super();
        this.areFileOperationsActivated = false;
        this.outputDirectory = null;
        this.exceptionsPrintWriter = null;
        this.filteredMoleculesPrintWriter = null;
        this.dataOutputPrintWriter = null;
    }
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Public test methods">
    
    //<editor-fold defaultstate="collapsed" desc="Tests involving databases">
    /**
     * Test for analyzing molecules in an SD file for all four different electron donation models supplied by the cdk: 
     * daylight, cdk, piBonds, cdkAllowingExocyclic and the aromaticity model cdkLegacy.
     * <p>
     * (Functional groups occurring multiple times in the same molecule are counted multiple times)
     *
     * @throws java.lang.Exception if initializeWithFileOperations() throws an exception or an unexpected exception occurs
     */
    @Test
    public void testElectronDonationDependency() throws Exception {
        this.analyzeElectronDonationDependency(ErtlFunctionalGroupsFinderEvaluationTest.SD_FILE_PATH, 
                ErtlFunctionalGroupsFinderEvaluationTest.ELECTRON_DONATION_TEST_IDENTIFIER, 
                true);
    }
    
    /**
     * Test for analyzing molecules in an SD file for all four different electron donation models supplied by the cdk: 
     * daylight, cdk, piBonds, cdkAllowingExocyclic and the aromaticity model cdkLegacy.
     * <p>
     * Difference to testElectronDonationDependency(): If the same functional group occurs multiple times in the same molecule
     * it is counted only once
     *
     * @throws java.lang.Exception if initializeWithFileOperations() throws an exception or an unexpected exception occurs
     */
    @Test
    public void testElectronDonationDependencyNoMultiples() throws Exception {
        this.analyzeElectronDonationDependency(ErtlFunctionalGroupsFinderEvaluationTest.SD_FILE_PATH, 
                ErtlFunctionalGroupsFinderEvaluationTest.ELECTRON_DONATION_NO_MULTIPLES_TEST_IDENTIFIER, 
                false);
    }
    
    /**
     * Test for analyzing molecules in an SD file for six different CycleFinder settings supplied by the cdk: all(), 
     * vertexShort(), relevant(), essential(), tripleShort() and cdkAromaticSet().
     * <p>
     * (Functional groups occurring multiple times in the same molecule are counted multiple times)
     *
     * @throws java.lang.Exception if initializeWithFileOperations() throws an exception or an unexpected exception occurs
     */
    @Test
    public void testCycleFinderDependency() throws Exception {
        this.initializeWithFileOperations(ErtlFunctionalGroupsFinderEvaluationTest.SD_FILE_PATH, 
                ErtlFunctionalGroupsFinderEvaluationTest.CYCLE_FINDER_TEST_IDENTIFIER);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("\nLoading file with path: " + ErtlFunctionalGroupsFinderEvaluationTest.SD_FILE_PATH);
        File tmpSDFile = new File(ErtlFunctionalGroupsFinderEvaluationTest.SD_FILE_PATH);
        int tmpRequiredNumberOfReaders = 6;
        IteratingSDFReader[] tmpReaders = new IteratingSDFReader[tmpRequiredNumberOfReaders];
        try {
            for (int i = 0; i < tmpRequiredNumberOfReaders; i++) {
                IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile),
                        DefaultChemObjectBuilder.getInstance(), true);
                tmpReaders[i] = tmpReader;
            }
        } catch (FileNotFoundException aFileNotFoundException) {
            System.out.println("\nSD file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
            return;
        }
        
        Aromaticity tmpDaylightModelAll = new Aromaticity(ElectronDonation.daylight(), Cycles.all());
        Aromaticity tmpDaylightModelVertexShort = new Aromaticity(ElectronDonation.daylight(), Cycles.vertexShort());
        Aromaticity tmpDaylightModelRelevant = new Aromaticity(ElectronDonation.daylight(), Cycles.relevant());
        Aromaticity tmpDaylightModelEssential = new Aromaticity(ElectronDonation.daylight(), Cycles.essential());
        Aromaticity tmpDaylightModelTripleShort = new Aromaticity(ElectronDonation.daylight(), Cycles.tripletShort());
        Aromaticity tmpDaylightModelCdkAromaticSet = new Aromaticity(ElectronDonation.daylight(), Cycles.cdkAromaticSet());
        
        boolean tmpAreMultiplesCounted = false;
        
        this.calculateAbsoluteFGFrequencies(tmpReaders[0], "all", tmpDaylightModelAll, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[1], "vertexShort", tmpDaylightModelVertexShort, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[2], "relevant", tmpDaylightModelRelevant, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[3], "essential", tmpDaylightModelEssential, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[4], "tripletShort", tmpDaylightModelTripleShort, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[5], "cdkAromaticSet", tmpDaylightModelCdkAromaticSet, tmpAreMultiplesCounted);
        
        System.out.println("\nAll analyses are done!");
        for (IteratingSDFReader tmpReader : tmpReaders) {
            tmpReader.close();
        }
        this.saveData();
        System.out.println("\nFinished!");
        System.out.println("\nNumber of occured exceptions: " + this.exceptionsCounter);
    }
    
    /**
     * Testing the ErtlFunctionalGroupsFinder.find() method's performance on the given SD file.
     * 
     * @throws java.lang.Exception if initializeWithFileOperations() throws an exception or the IteratingSDFReader 
     * can not be closed or an unexpectedException occurs
     */
    @Test
    public void testPerformance() throws Exception {
        this.initialize(true, "PerformanceTest");
        //First, check if the SD file is present and ignore test if it is not
        String tmpPathToSDFile = ErtlFunctionalGroupsFinderEvaluationTest.SD_FILE_PATH;
        System.out.println("\nLoading file with path: " + tmpPathToSDFile);
        File tmpSDFile = new File(tmpPathToSDFile);
        if (!tmpSDFile.canRead()) {
            System.out.println("\n\tUnable to find or read a file with path \"" + tmpPathToSDFile + "\".");
            System.out.println("\nTest is ignored.");
            Assume.assumeTrue(false);
            return;
        }
        IteratingSDFReader tmpReader;
        try {
            tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile),
                    DefaultChemObjectBuilder.getInstance(), true);
        } catch (FileNotFoundException aFileNotFoundException) {
            System.out.println("\nSD file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
            return;
        }
        List<IAtomContainer> tmpMoleculesList = new LinkedList<>();
        Aromaticity tmpCdkLegacyModel = new Aromaticity(ElectronDonation.daylight(), 
                Cycles.or(Cycles.all(), Cycles.vertexShort()));
        while (tmpReader.hasNext()) {
            try {
                IAtomContainer tmpMolecule = (IAtomContainer) tmpReader.next();
                tmpMolecule = this.applyFiltersAndPreprocessing(tmpMolecule);
                if (tmpMolecule.getProperty(ErtlFunctionalGroupsFinderEvaluationTest.MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY)) {
                    /*No logging required here, it is a simple performance test*/
                    continue;
                }
                tmpCdkLegacyModel.apply(tmpMolecule);
                tmpMoleculesList.add(tmpMolecule);
            } catch (Exception anException) {
                /*No logging required here, it is a simple performance test*/
            }
        }
        tmpReader.close();
        IAtomContainer[] tmpMoleculesArray = new IAtomContainer[tmpMoleculesList.size()];
        tmpMoleculesArray = tmpMoleculesList.toArray(tmpMoleculesArray);
        System.out.println("\nDone Loading database. Found and processed " + tmpMoleculesArray.length + " valid molecules.");
        long tmpStartTime = System.currentTimeMillis();
        for (IAtomContainer tmpMolecule : tmpMoleculesArray) {
            this.ertlFGFinderGenOn.find(tmpMolecule, false);
        }
        long tmpEndTime = System.currentTimeMillis();
        System.out.println("\nExtraction of functional groups from these molecules took " + (tmpEndTime - tmpStartTime) + " ms.\n");
    }
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Other tests">
    /**
     * Test for correct MoleculeHashGenerator settings/performance on some examples.
     *
     * @throws java.lang.Exception if initialize() throws an exception or a SMILES code can not be parsed into a molecule
     */
    @Test
    public void testMoleculeHashGeneratorSettings() throws Exception {
        this.initialize(false, "");
        SmilesParser tmpSmilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        
        /*Chebi70986, Chebi16238 and Chebi57692 all contain the same functional group with pseudo SMILES code
        "O=C1N=C(C(=NR)C(=O)N1R)N(R)R", but different hybridizations in the resulting atom containers. But their hash
        codes should be the same under the given settings. This is tested exemplary for many similar cases*/
        String[] tmpSmilesArray = {"OC[C@@H](O)[C@@H](O)[C@@H](O)CN1CC(CO)N=C2C(=O)NC(=O)N=C12",
            "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)OP(O)(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)n3cnc4c(N)ncnc34)c2cc1C",
            "Cc1cc2nc3c(nc(=O)[n-]c3=O)n(C[C@H](O)[C@H](O)[C@H](O)COP([O-])(=O)OP([O-])(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)n3cnc4c(N)ncnc34)c2cc1C"};
        List<Long> tmpHashCodesList = new LinkedList<>();
        for (String tmpSmilesCode : tmpSmilesArray) {
            IAtomContainer tmpParsedMolecule = tmpSmilesParser.parseSmiles(tmpSmilesCode);
            tmpParsedMolecule = this.applyFiltersAndPreprocessing(tmpParsedMolecule);
            Aromaticity.cdkLegacy().apply(tmpParsedMolecule);
            List<IAtomContainer> tmpFunctionalGroups = this.ertlFGFinderGenOn.find(tmpParsedMolecule);
            for (IAtomContainer tmpFunctionalGroup : tmpFunctionalGroups) {
                if (this.getPseudoSmilesCode(tmpFunctionalGroup).equals("O=C1N=C(C(=NR)C(=O)N1R)N(R)R")) {
                    tmpHashCodesList.add(this.molHashGenerator.generate(tmpFunctionalGroup));
                }
            }
        }
        for (Long tmpHashCode1 : tmpHashCodesList) {
            for (Long tmpHashCode2 : tmpHashCodesList) {
                Assert.assertEquals(tmpHashCode1.longValue(), tmpHashCode2.longValue());
            }
        }
        
        /*Functional groups like the tertiary amine or the hydroxyl group appear with aromatic and non-aromatic central
        atoms. These two cases should be discrimated by the MoleculeHashGenerator under the given settings*/
        String tmpTertiaryAmineSmiles = "*N(*)*";
        IAtomContainer tmpAromMol = tmpSmilesParser.parseSmiles(tmpTertiaryAmineSmiles);
        IAtomContainer tmpNonAromMol = tmpSmilesParser.parseSmiles(tmpTertiaryAmineSmiles);
        for (IAtom tmpAtom : tmpAromMol.atoms()) {
            if (tmpAtom.getSymbol().equals("N"))
                tmpAtom.setIsAromatic(true);
        }
        Assert.assertNotEquals(this.molHashGenerator.generate(tmpAromMol), this.molHashGenerator.generate(tmpNonAromMol));
        String tmpHydroxylGroupSmiles = "[H]O[C]";
        tmpAromMol = tmpSmilesParser.parseSmiles(tmpHydroxylGroupSmiles);
        tmpNonAromMol = tmpSmilesParser.parseSmiles(tmpHydroxylGroupSmiles);
        for (IAtom tmpAtom : tmpAromMol.atoms()) {
            if (tmpAtom.getSymbol().equals("C"))
                tmpAtom.setIsAromatic(true);
        }
        Assert.assertNotEquals(this.molHashGenerator.generate(tmpAromMol), this.molHashGenerator.generate(tmpNonAromMol));
        
        /*The following are examples of different (unique!) SMILES codes representing the same functional groups.
        They should be assigned the same hash code*/
        HashMap<String,String> tmpEquivalentSmilesMap = new HashMap<>(20);
        tmpEquivalentSmilesMap.put("*[N](*)=C(N(*)*)N(*)*", "*N(*)C(=[N](*)*)N(*)*");
        tmpEquivalentSmilesMap.put("*SC1=[N](*)[C]=[C]N1*", "*SC=1N(*)[C]=[C][N]1*");
        tmpEquivalentSmilesMap.put("*[N]1=[C][C]=[C]N1*", "*N1[C]=[C][C]=[N]1*");
        tmpEquivalentSmilesMap.put("*[N](*)=[C]N(*)*", "*N(*)[C]=[N](*)*");
        tmpEquivalentSmilesMap.put("*N(*)[C]=[C][C]=[C][C]=[C][C]=[C][C]=[N](*)*", "*[N](*)=[C][C]=[C][C]=[C][C]=[C][C]=[C]N(*)*");
        tmpEquivalentSmilesMap.put("*[N](*)=C(N(*)*)N(*)P(=O)(O[H])O[H]", "*N(*)C(=[N](*)*)N(*)P(=O)(O[H])O[H]");
        tmpEquivalentSmilesMap.put("[O]I(=O)=O", "O=I(=O)[O]");
        tmpEquivalentSmilesMap.put("[O]Br(=O)=O", "O=Br(=O)[O]");
        tmpEquivalentSmilesMap.put("[O]Cl(=O)(=O)=O", "O=Cl(=O)(=O)[O]");
        tmpEquivalentSmilesMap.put("[C]=[C][C]=[C]C#C[C]=[C]C#[C]", "[C]#C[C]=[C]C#C[C]=[C][C]=[C]");
        tmpEquivalentSmilesMap.put("*N1[C]=[C][C]=[N]1*", "*[N]1=[C][C]=[C]N1*");
        tmpEquivalentSmilesMap.put("O=C(*)O*", "*OC(*)=O");
        for (String tmpKeySmiles : tmpEquivalentSmilesMap.keySet()) {
            IAtomContainer tmpKeyMol = tmpSmilesParser.parseSmiles(tmpKeySmiles);
            IAtomContainer tmpValueMol = tmpSmilesParser.parseSmiles(tmpEquivalentSmilesMap.get(tmpKeySmiles));
            Assert.assertEquals(this.molHashGenerator.generate(tmpKeyMol), this.molHashGenerator.generate(tmpValueMol));
        }
    }
    
    /**
     * Test for correct preprocessing (neutralization of charges and selection of biggest fragment).
     * 
     * @throws Exception if initialize() throws an exception or a SMILES code can not be parsed into a molecule
     */
    @Test
    public void testPreprocessing() throws Exception {
        this.initialize(false, "");
        String tmpSmiles = "CC[O-].C";
    	SmilesParser tmpSmilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMol = tmpSmilesParser.parseSmiles(tmpSmiles);
        tmpMol = this.applyFiltersAndPreprocessing(tmpMol);
	SmilesGenerator tmpGenerator = SmilesGenerator.unique();
	Assert.assertEquals("OCC", tmpGenerator.create(tmpMol));
    }
    //</editor-fold>
    
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Analyzes molecules in an SD file for all four different electron donation models supplied by the cdk: 
     * daylight, cdk, piBonds, cdkAllowingExocyclic and the aromaticity model cdkLegacy.
     * 
     * @param anSDFilePath absolute path of the SD file to analyze
     * @param aTestIdentifier a folder with this name will be created in the output directory and it will be added to 
     * the output and log files' names for association of test and files; may be null or empty
     * @param anAreMultiplesCounted if false, functional groups that occur multiple times in the same molecule will 
     * only be counted once
     * @throws java.lang.Exception if initializeWithFileOperations() throws an exception or an unexpected exception occurs
     */
    private void analyzeElectronDonationDependency(String anSDFilePath, 
            String aTestIdentifier, 
            boolean anAreMultiplesCounted) throws Exception {
        this.initializeWithFileOperations(anSDFilePath, aTestIdentifier);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("\nLoading file with path: " + anSDFilePath);
        File tmpSDFile = new File(anSDFilePath);
        int tmpRequiredNumberOfReaders = 5;
        IteratingSDFReader[] tmpReaders = new IteratingSDFReader[tmpRequiredNumberOfReaders];
        try {
            for (int i = 0; i < tmpRequiredNumberOfReaders; i++) {
                IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile),
                        DefaultChemObjectBuilder.getInstance(), true);
                tmpReaders[i] = tmpReader;
            }
        } catch (FileNotFoundException aFileNotFoundException) {
            System.out.println("\nSD file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
            return;
        }
        //If the 'all' CycleFinder produces an Intractable exception the 'vertexShort' CycleFinder is used
        CycleFinder tmpCycleFinder = Cycles.or(Cycles.all(), Cycles.vertexShort());
        
        Aromaticity tmpDaylightModel = new Aromaticity(ElectronDonation.daylight(), tmpCycleFinder);
        Aromaticity tmpCdkModel = new Aromaticity(ElectronDonation.cdk(), tmpCycleFinder);
        Aromaticity tmpPiBondsModel = new Aromaticity(ElectronDonation.piBonds(), tmpCycleFinder);
        Aromaticity tmpCdkAllowingExocyclicModel = new Aromaticity(ElectronDonation.cdkAllowingExocyclic(), tmpCycleFinder);
        Aromaticity tmpCDKLegacyModel = Aromaticity.cdkLegacy();
        
        this.calculateAbsoluteFGFrequencies(tmpReaders[0], 
                ErtlFunctionalGroupsFinderEvaluationTest.DAYLIGHT_MODEL_SETTINGS_KEY, tmpDaylightModel, anAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[1], 
                ErtlFunctionalGroupsFinderEvaluationTest.CDK_MODEL_SETTINGS_KEY, tmpCdkModel, anAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[2], 
                ErtlFunctionalGroupsFinderEvaluationTest.PIBONDS_MODEL_SETTINGS_KEY, tmpPiBondsModel, anAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[3], 
                ErtlFunctionalGroupsFinderEvaluationTest.CDK_EXOCYCLIC_MODEL_SETTINGS_KEY, tmpCdkAllowingExocyclicModel, anAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[4], 
                ErtlFunctionalGroupsFinderEvaluationTest.CDK_LEGACY_MODEL_SETTINGS_KEY, tmpCDKLegacyModel, anAreMultiplesCounted);
        
        System.out.println("\nAll analyses are done!");
        for (IteratingSDFReader tmpReader : tmpReaders) {
            tmpReader.close();
        }
        this.saveData();
        System.out.println("\nFinished!");
        System.out.println("\nNumber of occured exceptions: " + this.exceptionsCounter);
    }
    
    /**
     * Initializes all class variables except the working directory and the PrintWriter instances. This method should be 
     * called directly when a test does not require any of the specific file operations like logging or reading an SD file.
     * 
     * @param aShouldPrintHeader true, if this method is called directly (not by initializeWithFileOperations()) and 
     * should print on the console that a new test was started
     * @param aTestIdentifier if aShouldPrintHeader is true this ID will be printed to the console
     */
    private void initialize(boolean aShouldPrintHeader, String aTestIdentifier) {
        if (aShouldPrintHeader) {
            System.out.println("\n#########################################################################\n");
            System.out.println("Starting new test, identifier: " + aTestIdentifier);
            System.out.println("\nInitializing class variables...");
        }
        this.smilesGenerator = new SmilesGenerator(ErtlFunctionalGroupsFinderEvaluationTest.SMILES_GENERATOR_OUTPUT_MODE);
        this.molHashGenerator = new HashGeneratorMaker()
                .depth(8)
                .elemental()
                //following line is used instead of .orbital() because the atom hybridizations take more information into 
                //account than the bond order sum but that is not required here
                //Note: This works here because the ErtlFunctionalGroupsFinder extracts the relevant atoms and bonds only
                //resulting in incomplete valences that can be used here in this way
                .encode(BasicAtomEncoder.BOND_ORDER_SUM)
                .encode(CustomAtomEncoder.AROMATICITY) //See enum CustomAtomEncoder below
                .molecular();
        this.ertlFGFinderGenOff = new ErtlFunctionalGroupsFinder(Mode.NO_GENERALIZATION);
        this.ertlFGFinderGenOn = new ErtlFunctionalGroupsFinder(Mode.DEFAULT);
        this.masterHashMap = new HashMap(ErtlFunctionalGroupsFinderEvaluationTest.MASTER_HASHMAP_INITIAL_CAPACITY, 
                ErtlFunctionalGroupsFinderEvaluationTest.MASTER_HASHMAP_LOAD_FACTOR);
        this.settingsKeysList = new LinkedList<>();
        this.exceptionsCounter = 0;
        this.areFilteredMoleculesLogged = false;
        String[] tmpMetalNumbersStrings = ErtlFunctionalGroupsFinderEvaluationTest.NON_METALLIC_ATOMIC_NUMBERS.split(",");
        Integer[] tmpMetalNumbersInt = new Integer[tmpMetalNumbersStrings.length];
        for (int i = 0; i < tmpMetalNumbersStrings.length; i++) {
            tmpMetalNumbersInt[i] = Integer.parseInt(tmpMetalNumbersStrings[i]);
        }
        this.nonMetallicAtomicNumbersSet = new HashSet(Arrays.asList(tmpMetalNumbersInt));
        this.pseudoSmilesAromaticElementToPlaceholderElementMap = new HashMap<>(10, 1);
        this.pseudoSmilesAromaticElementToPlaceholderElementMap.put("C", "Ce");
        this.pseudoSmilesAromaticElementToPlaceholderElementMap.put("N", "Nd");
        this.pseudoSmilesAromaticElementToPlaceholderElementMap.put("S", "Sm");
        this.pseudoSmilesAromaticElementToPlaceholderElementMap.put("O", "Os");
        this.pseudoSmilesAromaticElementToPlaceholderElementMap.put("Se", "Sc");
        this.pseudoSmilesAromaticElementToPlaceholderElementMap.put("P", "Pm");
        this.pseudoSmilesAromaticElementToPlaceholderElementMap.put("R", "Es");
        this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap = new HashMap<>(10, 1);
        this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.put("Es", 
                ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_R_ATOM);
        this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.put("Pm", 
                ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_AROMATIC_PHOSPHOR);
        this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.put("Sc", 
                ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_AROMATIC_SELENIUM);
        this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.put("Os", 
                ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_AROMATIC_OXYGEN);
        this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.put("Sm",
                ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_AROMATIC_SULPHUR);
        this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.put("Nd", 
                ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_AROMATIC_NITROGEN);
        this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.put("Ce", 
                ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_AROMATIC_CARBON);
        this.isTestAbleToRun = true;
        if (aShouldPrintHeader) {
            System.out.println("\n\tDone initializing class variables.\n");
        }
    }
    /**
     * Initializes all class variables and determines the output directory.
     * 
     * @param anSDFilePath absolute path of the SD file to analyze for a quick pre-check if it is present and the test 
     * is therefore meant to run; may be empty but not null
     * @param aTestIdentifier a folder with this name will be created in the output directory and it will be added to 
     * the output and log files' names for association of test and files; may be null or empty
     * @throws java.lang.Exception if one the FileWriter instances can not be instantiated, more than 
     * Integer.MAX-VALUE tests are to be run this minute (error in the naming of output files), aPathOfSDFile is null or 
     * an unexpected exception occurs.
     */
    private void initializeWithFileOperations(String anSDFilePath, String aTestIdentifier) throws Exception {
        System.out.println("\n#########################################################################\n");
        System.out.println("Starting new test, identifier: " + aTestIdentifier);
        System.out.println("\nInitializing class variables...");
        this.isTestAbleToRun = true;
        //First, check if the SD file is present and ignore test if it is not
        File tmpSDFile = new File(anSDFilePath);
        if (!tmpSDFile.canRead() || tmpSDFile.getAbsoluteFile().getParent() == null) {
            System.out.println("\n\tUnable to find or read a file with path \"" + anSDFilePath + "\" or to get its parent directory.");
            System.out.println("\nTest is ignored.");
            this.isTestAbleToRun = false;
            Assume.assumeTrue(false);
            return;
        }
        //Determine the output directory
        String tmpOutputRootDirectory = tmpSDFile.getAbsoluteFile().getParent() + File.separator;
        this.outputDirectory = tmpOutputRootDirectory 
                + ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FOLDER_FROM_WORKING_DIRECTORY 
                + File.separator 
                + aTestIdentifier;
        File tmpOutputDirectoryFile = new File(this.outputDirectory);
        if (!tmpOutputDirectoryFile.exists()) {
            tmpOutputDirectoryFile.mkdirs();
        }
        System.out.println("\n\tOutput directory: " + this.outputDirectory);
        //Create a time stamp for output and log files
        LocalDateTime tmpDateTime = LocalDateTime.now();
        String tmpDateTimeAddition = tmpDateTime.format(DateTimeFormatter.ofPattern(
                ErtlFunctionalGroupsFinderEvaluationTest.DATE_TIME_FORMAT_PATTERN));
        //Set up exceptions log file
        File tmpExceptionsLogFile = new File(this.outputDirectory + File.separator + aTestIdentifier 
                + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                + ErtlFunctionalGroupsFinderEvaluationTest.EXCEPTIONS_LOG_FILE_NAME 
                + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                + tmpDateTimeAddition 
                + ErtlFunctionalGroupsFinderEvaluationTest.EXCEPTIONS_LOG_FILE_TYPE);
        int tmpFilesInThisMinuteCounter = 1;
        boolean tmpNumberAddedToFileName = false;
        //No pre-existing file is overwritten
        if (tmpExceptionsLogFile.exists()) {
            tmpNumberAddedToFileName = true;
            while (tmpFilesInThisMinuteCounter <= Integer.MAX_VALUE) {
                tmpExceptionsLogFile = new File(this.outputDirectory + File.separator + aTestIdentifier 
                        + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                        + ErtlFunctionalGroupsFinderEvaluationTest.EXCEPTIONS_LOG_FILE_NAME 
                        + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                        + tmpDateTimeAddition 
                        + "(" + tmpFilesInThisMinuteCounter + ")" 
                        + ErtlFunctionalGroupsFinderEvaluationTest.EXCEPTIONS_LOG_FILE_TYPE);
                if (!tmpExceptionsLogFile.exists()) {
                    break;
                }
                if (tmpFilesInThisMinuteCounter == Integer.MAX_VALUE) {
                    throw new Exception("More than [Integer.MAX-VALUE] tests are to be run this minute. "
                            + "This test class is not configurated for that.");
                }
                tmpFilesInThisMinuteCounter++;
            }
        }
        FileWriter tmpExceptionsLogFileWriter = new FileWriter(tmpExceptionsLogFile);
        this.exceptionsPrintWriter = new PrintWriter(tmpExceptionsLogFileWriter);
        this.exceptionsPrintWriter.println(ErtlFunctionalGroupsFinderEvaluationTest.EXCEPTIONS_LOG_FILE_HEADER);
        this.exceptionsPrintWriter.flush();
        System.out.println("\tExceptions will be written to: " + tmpExceptionsLogFile.getName());
        //Set up filtered molecules log file
        File tmpFilteredMoleculesFile;
        if (tmpNumberAddedToFileName) {
            tmpFilteredMoleculesFile = new File(this.outputDirectory+ File.separator + aTestIdentifier 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILTERED_MOLECULES_FILE_NAME 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                    + tmpDateTimeAddition + "(" + tmpFilesInThisMinuteCounter + ")" 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILTERED_MOLECULES_FILE_TYPE);
        } else {
            tmpFilteredMoleculesFile = new File(this.outputDirectory+ File.separator + aTestIdentifier 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILTERED_MOLECULES_FILE_NAME 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                    + tmpDateTimeAddition 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILTERED_MOLECULES_FILE_TYPE);
        }
        FileWriter tmpFilteredMoleculesFileWriter = new FileWriter(tmpFilteredMoleculesFile);
        this.filteredMoleculesPrintWriter = new PrintWriter(tmpFilteredMoleculesFileWriter);
        this.filteredMoleculesPrintWriter.println(ErtlFunctionalGroupsFinderEvaluationTest.FILTERED_MOLECULES_FILE_HEADER);
        this.filteredMoleculesPrintWriter.flush();
        System.out.println("\tFiltered molecules will be written to: " + tmpFilteredMoleculesFile.getName());
        //Set up output file
        File tmpOutputFile;
        if (tmpNumberAddedToFileName) {
            tmpOutputFile = new File(this.outputDirectory+ File.separator + aTestIdentifier 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                    + ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_NAME 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                    + tmpDateTimeAddition 
                    + "(" + tmpFilesInThisMinuteCounter + ")" 
                    + ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_TYPE);
        } else {
            tmpOutputFile = new File(this.outputDirectory+ File.separator + aTestIdentifier 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR
                    + ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_NAME 
                    + ErtlFunctionalGroupsFinderEvaluationTest.FILE_NAME_ADDITION_SEPERATOR 
                    + tmpDateTimeAddition 
                    + ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_TYPE);
        }
        FileWriter tmpOutputFileWriter = new FileWriter(tmpOutputFile);
        this.dataOutputPrintWriter = new PrintWriter(tmpOutputFileWriter);
        System.out.println("\tThe absolute functional groups frequencies will be written to: " + tmpOutputFile.getName());
        this.areFileOperationsActivated = true;
        this.initialize(false, aTestIdentifier);
    }
    
    /**
     * Does one iteration of loading molecules from an IteratingSDFReader, applying the given aromaticity model, 
     * extracting their functional groups (with and without generalization) and adding the results to the master HashMap.
     * Exceptions caused by read-in molecules are caught and logged.
     * 
     * @param aReader to load the molecules to be screened from an SD file
     * @param aSettingsKey resulting functional groups will be added to the inner maps of the master HasMap under this key
     * @param anAromaticity to apply to the molecules
     * @param anAreMultiplesCounted if false, functional groups that occur multiple times in the same molecule will 
     * only be counted once
     */
    private void calculateAbsoluteFGFrequencies(
            IteratingSDFReader aReader, 
            String aSettingsKey, 
            Aromaticity anAromaticity, 
            boolean anAreMultiplesCounted) {
        
        System.out.println("\nAnalyzing database using specified settings: " + aSettingsKey);
        //<editor-fold defaultstate="collapsed" desc="Counter definitions">
        int tmpMoleculesCounter = 0; //total number of molecules successfully loaded from the Sd file
        int tmpChargeCounter = 0; //number of molecules that contain one or more charged atoms
        int tmpUnconnectedCounter = 0; //number of structures with unconnected substructures
        int tmpFilteredMoleculesCounter = 0; //filtered molecules, see below
        int tmpMetallicCounter = 0; //number of molecules that were filtered because they contain one or more metal, metalloid or "R" atoms
        int tmpNoAtomOrBondCounter = 0; //number of molecules that were filtered because they do not contain an atom or a bond
        int tmpSkippedMoleculesCounter = 0; //molecules that caused exceptions and are therefore skipped
        int tmpValidMoleculesCounter = 0; //molecules that were successfully passed to the find() method without causing an exception
        int tmpNoFunctionalGroupsCounter = 0; //number of molecules that contain no functional group
        //</editor-fold>
        while (aReader.hasNext()) {
            List<IAtomContainer> tmpFunctionalGroups;
            List<IAtomContainer> tmpFunctionalGroupsGeneralized;
            IAtomContainer tmpMolecule = (IAtomContainer) aReader.next();
            tmpMoleculesCounter++;
            if (tmpMolecule == null) {
                break;
            }
            String tmpSettingsKeyForLogging = aSettingsKey;
            IAtomContainer tmpOriginalMolecule = null;
            try {
                /*Note: A molecule can be supposed to be filtered for more than one of the named reasons but only the first 
                  tested reason will be named as cause*/
                String tmpCauseForFiltering;
                //Remove s-groups before cloning, they cause trouble in atomContainer.clone(); see cdk api CDKConstants.CTAB_SGROUPS
                tmpMolecule.removeProperty(CDKConstants.CTAB_SGROUPS);
                //Produce a clone of the original molecule for logging
                tmpOriginalMolecule = tmpMolecule.clone();
                //Apply filters and preprocessing
                tmpMolecule = this.applyFiltersAndPreprocessing(tmpMolecule);
                //Filter molecule if necessary
                if (tmpMolecule.getProperty(ErtlFunctionalGroupsFinderEvaluationTest.MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY)) {
                    tmpCauseForFiltering = tmpMolecule.getProperty(
                            ErtlFunctionalGroupsFinderEvaluationTest.CAUSE_FOR_FILTERING_PROPERTY_KEY);
                    if (tmpCauseForFiltering.equals(ErtlFunctionalGroupsFinderEvaluationTest.ATOM_OR_BOND_COUNT_ZERO)) {
                        tmpNoAtomOrBondCounter++;
                    } else if (tmpCauseForFiltering.equals(ErtlFunctionalGroupsFinderEvaluationTest.FORBIDDEN_ATOMIC_NUMBER)) {
                        tmpMetallicCounter++;
                    }
                    tmpFilteredMoleculesCounter++;
                    if (!this.areFilteredMoleculesLogged && this.areFileOperationsActivated) {
                        this.logFilteredMolecule(tmpMolecule, tmpFilteredMoleculesCounter, tmpCauseForFiltering);
                    }
                    continue;
                }
                //Add according additions to the settings key for logging if the molecule was preprocessed
                if (tmpMolecule.getProperty(ErtlFunctionalGroupsFinderEvaluationTest.BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY)) {
                    tmpUnconnectedCounter++;
                    tmpSettingsKeyForLogging += ErtlFunctionalGroupsFinderEvaluationTest.FRAGMENT_SELECTED_SETTINGS_KEY_ADDITION;
                }
                if (tmpMolecule.getProperty(ErtlFunctionalGroupsFinderEvaluationTest.CHARGES_NEUTRALIZED_PROPERTY_KEY)) {
                    tmpChargeCounter++;
                    tmpSettingsKeyForLogging += ErtlFunctionalGroupsFinderEvaluationTest.NEUTRALIZED_SETTINGS_KEY_ADDITION;
                }
                //Application of aromaticity model
                anAromaticity.apply(tmpMolecule);
                //Now the analysis of functional groups
                //In the first invokation of find() the atom container is cloned because it will be reused for the second invokation of find()
                tmpFunctionalGroups = this.ertlFGFinderGenOff.find(tmpMolecule, true);
                //Do the extraction again with activated generalization
                tmpSettingsKeyForLogging += ErtlFunctionalGroupsFinderEvaluationTest.GENERALIZATION_SETTINGS_KEY_ADDITION;
                //In this invokation of find() the atom container is not cloned
                tmpFunctionalGroupsGeneralized = this.ertlFGFinderGenOn.find(tmpMolecule, false);
                tmpValidMoleculesCounter++;
            } catch (Exception anException) {
                tmpSkippedMoleculesCounter++;
                if (this.areFileOperationsActivated) {
                    if (tmpOriginalMolecule != null) {
                        this.logException(anException, tmpSettingsKeyForLogging, tmpOriginalMolecule);
                    } else {
                        this.logException(anException, tmpSettingsKeyForLogging, tmpMolecule);
                    }
                }
                continue;
            }
            if (!(tmpFunctionalGroups == null || tmpFunctionalGroups.isEmpty())) {
                //If a molecule does not have FGs without generalization it does not have FGs with generalization either
                this.enterFunctionalGroupsIntoMasterMap(tmpFunctionalGroups, 
                        aSettingsKey, 
                        anAreMultiplesCounted, 
                        tmpOriginalMolecule);
                this.enterFunctionalGroupsIntoMasterMap(tmpFunctionalGroupsGeneralized, 
                        aSettingsKey + ErtlFunctionalGroupsFinderEvaluationTest.GENERALIZATION_SETTINGS_KEY_ADDITION,
                        anAreMultiplesCounted, 
                        tmpOriginalMolecule);
            } else {
                tmpNoFunctionalGroupsCounter++;
            }
        }
        try {
            aReader.close();
        } catch (IOException anIOException) { }
        //Since the filters remain the same in every iteration filtered molecules must be logged only once
        //(assuming that only one SD file is analyzed in a test)
        if (!this.areFilteredMoleculesLogged) {
            this.areFilteredMoleculesLogged = true;
        }
        System.out.println("Analysis done");
        System.out.println(tmpMoleculesCounter + " molcules were screened.");
        System.out.println(tmpChargeCounter + " molecules contained one or more charged atom.");
        System.out.println(tmpUnconnectedCounter + " molecules had two or more unconnected structures.");
        System.out.println(tmpFilteredMoleculesCounter + " molecules were filtered.");
        System.out.println("\t" + tmpNoAtomOrBondCounter + " molecules were filtered because their atom or bond count is 0.");
        System.out.println("\t" + tmpMetallicCounter + " molecules were filtered because they contained one or more "
                + "metal, metalloid or \"R\" atom.");
        System.out.println(tmpSkippedMoleculesCounter + " molecules were skipped due to exceptions.");
        System.out.println(tmpValidMoleculesCounter + " molecules were valid.");
        System.out.println(tmpNoFunctionalGroupsCounter + " molecules contained no functional groups.");
        System.out.println(this.masterHashMap.size() + " different functional groups were detected so far.");
        
    }
    
    /**
     * Combines all filtering and preprocessing steps. Molecules will be filtered if they contain elements with a not allowed atomic 
     * number (metal, metalloid or 'R' atoms) or if their atom or bond count is zero. If the molecule should be filtered
     * the IAtomContainer property MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY will be set to true and the 
     * CAUSE_FOR_FILTERING_PROPERTY_KEY will give the cause for the filtering.
     * <p>
     * The preprocessing consists of neutralizing any charges in the molecule and selecting the biggest fragment for
     * further processing if the molecule consists of one or more unconnected structures. If any of these cases apply
     * the IAtomContainer properties CHARGES_NEUTRALIZED_PROPERTY_KEY and BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY will be
     * set accordingly.
     * 
     * @param aMolecule the molecule to be processed
     * @return the processed molecule
     * @throws CDKException if AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms() or neutralizeCharges() 
     * throws a CDKException
     */
    private IAtomContainer applyFiltersAndPreprocessing(IAtomContainer aMolecule) throws CDKException {
        aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY, false);
        aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.CAUSE_FOR_FILTERING_PROPERTY_KEY, "");
        aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.CHARGES_NEUTRALIZED_PROPERTY_KEY, false);
        aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY, false);
        if (this.isAtomOrBondCountZero(aMolecule)) {
            aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY, true);
            aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.CAUSE_FOR_FILTERING_PROPERTY_KEY, 
                    ErtlFunctionalGroupsFinderEvaluationTest.ATOM_OR_BOND_COUNT_ZERO);
            return aMolecule;
        }
        /*Remove s-groups before cloning (will be done by the ErtlFunctionalGroupsFinder.find() method), 
        they cause trouble in atomContainer.clone(); see cdk api CDKConstants.CTAB_SGROUPS.*/
        aMolecule.removeProperty(CDKConstants.CTAB_SGROUPS);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
        //Preprocessing: From structures containing two or more unconnected structures (e.g. ions) 
        //choose the largest structure for analysis
        if (this.isStructureUnconnected(aMolecule)) {
            aMolecule = this.selectBiggestFragment(aMolecule);
            aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY, true);
        }
        //Filter molecules containing metals, metalloids or "R" atoms
        if (this.areMetallicOrMetalloidAtomsInMolecule(aMolecule)) {
            aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY, true);
            aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.CAUSE_FOR_FILTERING_PROPERTY_KEY, 
                    ErtlFunctionalGroupsFinderEvaluationTest.FORBIDDEN_ATOMIC_NUMBER);
            return aMolecule;
        }
        //Preprocessing: Neutralize charges
        if (this.isMoleculeCharged(aMolecule)) {
            aMolecule = this.neutralizeCharges(aMolecule);
            aMolecule.setProperty(ErtlFunctionalGroupsFinderEvaluationTest.CHARGES_NEUTRALIZED_PROPERTY_KEY, true);
        }
        return aMolecule;
    }
    
    /**
     * Returns true, if the atom or bond count of the molecule is zero. This is a cause for filtering the molecule.
     * 
     * @param aMolecule the molecule to be tested
     * @return true, if the atom or bond count of the molecule is zero
     */
    private boolean isAtomOrBondCountZero(IAtomContainer aMolecule) {
        return aMolecule.getAtomCount() == 0 || aMolecule.getBondCount() == 0;
    }
    
    /**
     * Returns true, if a not allowed atomic number is detected in the molecule. This is a cause for filtering the molecule.
     * 
     * @param aMolecule the molecule to be tested
     * @return true, if the molecule contains a not allowed element
     */
    private boolean areMetallicOrMetalloidAtomsInMolecule(IAtomContainer aMolecule) {
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (!this.nonMetallicAtomicNumbersSet.contains(tmpAtom.getAtomicNumber())) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Returns true, if the molecule consists of two or more unconnected structures.
     * 
     * @param aMolecule the molecule to be tested
     * @return true, if the molecule consists of two or more unconnected structures
     */
    private boolean isStructureUnconnected(IAtomContainer aMolecule) {
        return (!ConnectivityChecker.isConnected(aMolecule));
    }
    
    /**
     * Returns the biggest unconnected fragment/structure of the given molecule. To pre-check if the molecule has 
     * two or more unconnected structures use isStructureConnected(). All set properties of aMolecule will be copied to 
     * the returned IAtomContainer.
     * 
     * @param aMolecule the molecule whose biggest fragment should be found
     * @return the biggest unconnected fragment/structure of the given molecule
     */
    private IAtomContainer selectBiggestFragment(IAtomContainer aMolecule) {
        IAtomContainerSet tmpFragmentsSet = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        IAtomContainer tmpBiggestFragment = null;
        for (IAtomContainer tmpFragment : tmpFragmentsSet.atomContainers()) {
            if (tmpBiggestFragment == null || tmpBiggestFragment.getAtomCount() < tmpFragment.getAtomCount()) {
                tmpBiggestFragment = tmpFragment;
            }
        }
        tmpBiggestFragment.setProperties(aMolecule.getProperties());
        return tmpBiggestFragment;
    }
    
    /**
     * Returns true, if the molecule contains charged atoms.
     * 
     * @param aMolecule the molecule to be tested
     * @return true, if the molecule contains charged atoms
     */
    private boolean isMoleculeCharged(IAtomContainer aMolecule) {
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (tmpAtom.getFormalCharge() != 0) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Neutralizes all non-zero charges in the given molecule. To pre-check if the molecule has charged atoms use 
     * isMoleculeCharged().
     * 
     * @param aMolecule the molecule to be neutralized
     * @return the same IAtomContainer instance as aMolecule but with neutralized charges
     * @throws CDKException if CDKAtomTypeMatcher.findMatchingAtomType() or CDKHydrogenAdder.addImplicitHydrogens 
     * throws a CDKException
     */
    private IAtomContainer neutralizeCharges(IAtomContainer aMolecule) throws CDKException {
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (tmpAtom.getFormalCharge() != 0) {
                tmpAtom.setFormalCharge(0);
                CDKHydrogenAdder tmpHAdder = CDKHydrogenAdder.getInstance(aMolecule.getBuilder());
                CDKAtomTypeMatcher tmpMatcher = CDKAtomTypeMatcher.getInstance(aMolecule.getBuilder());
                IAtomType tmpMatchedType = tmpMatcher.findMatchingAtomType(aMolecule, tmpAtom);
                if (tmpMatchedType != null) {
                    AtomTypeManipulator.configure(tmpAtom, tmpMatchedType);
                }
                tmpHAdder.addImplicitHydrogens(aMolecule, tmpAtom);
            }
        }
        return aMolecule;
    }
    
    /**
     * Inserts a list of IAtomContainers (the functional groups of one molecule) into the master HashMap. If the 
     * functional group is already inserted (with this settings key) its frequency for the given settings key is raised 
     * by one or else a new inner HashMap will be created for it.
     * 
     * @param aFunctionalGroupsList the functional groups of one molecule to be inserted
     * @param aSettingsKey will be the key of the inner HashMap inside the master HashMap
     * @param anAreMultiplesCounted if false, functional groups that occur multiple times in aFunctionalGroupsList will 
     * only be entered once into the master Hashmap
     * @param anFGContainingMolecule the molecule from which the functional groups originated; will be added to the 
     * master Hashmap
     */
    private void enterFunctionalGroupsIntoMasterMap(
            List<IAtomContainer> aFunctionalGroupsList, 
            String aSettingsKey, 
            boolean anAreMultiplesCounted,
            IAtomContainer anFGContainingMolecule) {
        
        if (!this.settingsKeysList.contains(aSettingsKey)) {
            this.settingsKeysList.add(aSettingsKey);
        }
        List<Long> tmpAlreadyEnteredFGsForThisMol = new LinkedList<>();
        for (IAtomContainer tmpFunctionalGroup : aFunctionalGroupsList) {
            long tmpHashCode = this.molHashGenerator.generate(tmpFunctionalGroup);
            //Case: Multiples are counted only once and the functional group was already entered for this molecule
            if (!anAreMultiplesCounted && tmpAlreadyEnteredFGsForThisMol.contains(tmpHashCode)) {
                continue;
            }
            //Case: functional group is already in the master HashMap
            if (this.masterHashMap.containsKey(tmpHashCode)) {
                HashMap<String, Object> tmpInnerMap = (HashMap)this.masterHashMap.get(tmpHashCode);
                //And a key-value pair for this settings key is already present too -> raise frequency by one
                if (tmpInnerMap.containsKey(aSettingsKey)) {
                    int tmpFrequency = (int)tmpInnerMap.get(aSettingsKey);
                    tmpFrequency++;
                    tmpInnerMap.put(aSettingsKey, tmpFrequency);
                //there is no key-value pair for this settings key in the inner HashMap -> create one
                } else {
                    tmpInnerMap.put(aSettingsKey, 1);
                }
            //The functional group did not occur before -> create a new inner HashMap for this molecule
            } else {
                HashMap<String,Object> tmpNewInnerMap = new HashMap(
                        ErtlFunctionalGroupsFinderEvaluationTest.INNER_HASHMAPS_INITIAL_CAPACITY);
                tmpNewInnerMap.put(ErtlFunctionalGroupsFinderEvaluationTest.MOLECULE_OF_ORIGIN_KEY, anFGContainingMolecule);
                tmpNewInnerMap.put(aSettingsKey, 1);
                String tmpSmilesCode;
                String tmpPseudoSmilesCode;
                try {
                    //Creation of unique SMILES code
                    tmpSmilesCode = this.smilesGenerator.create(tmpFunctionalGroup);
                    //Creation of pseudo SMILES code
                    tmpPseudoSmilesCode = this.getPseudoSmilesCode(tmpFunctionalGroup);
                } catch (CDKException | NullPointerException | CloneNotSupportedException anException) {
                    if (this.areFileOperationsActivated) {
                        this.logException(anException, aSettingsKey + "Creating SMILES code", tmpFunctionalGroup);
                    }
                    tmpSmilesCode = ErtlFunctionalGroupsFinderEvaluationTest.SMILES_CODE_PLACEHOLDER;
                    tmpPseudoSmilesCode = ErtlFunctionalGroupsFinderEvaluationTest.SMILES_CODE_PLACEHOLDER;
                }
                tmpNewInnerMap.put(ErtlFunctionalGroupsFinderEvaluationTest.SMILES_CODE_KEY, tmpSmilesCode);
                tmpNewInnerMap.put(ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_CODE_KEY, tmpPseudoSmilesCode);
                this.masterHashMap.put(tmpHashCode, tmpNewInnerMap);
            }
            tmpAlreadyEnteredFGsForThisMol.add(tmpHashCode);
        }   
    }
    
    /**
     * Writes all frequency data with the respective hash code, SMILES code, pseudo SMILES code and the ChEBI or ChEMBL 
     * id or CDK title of an exemplary molecule that contains this functional group for all functional groups in the 
     * master HashMap to the output file.
     * <p>
     * Note: The IAtomContainer object stored in the master HashMap's inner maps are cloned in this method for pseudo 
     * SMILES creation. And all PrintWriter instances will be closed.
     */
    private void saveData() {
        if(!this.areFileOperationsActivated) {
            System.out.println("\nFile operations are not activated, invokation of saveData() is therefore not possible.");
            return;
        }
        System.out.println("\nWriting to file...");
        //Writing the output file's header
        String tmpFileHeader = ErtlFunctionalGroupsFinderEvaluationTest.HASH_CODE_KEY 
                + ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR 
                + ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_CODE_KEY 
                + ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR 
                + ErtlFunctionalGroupsFinderEvaluationTest.SMILES_CODE_KEY;
        for (String tmpSettingsKey : this.settingsKeysList) {
            tmpFileHeader += ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR + tmpSettingsKey;
        }
        tmpFileHeader += ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR 
                + ErtlFunctionalGroupsFinderEvaluationTest.MOLECULE_OF_ORIGIN_KEY;
        this.dataOutputPrintWriter.println(tmpFileHeader);
        this.dataOutputPrintWriter.flush();
        Iterator tmpFunctionalGroupsIterator = this.masterHashMap.keySet().iterator();
        //Iteration for all molecules in the master HashMap
        while (tmpFunctionalGroupsIterator.hasNext()) {
            long tmpHashCode = (long)tmpFunctionalGroupsIterator.next();
            HashMap tmpInnerMap = (HashMap)this.masterHashMap.get(tmpHashCode);
            String tmpSmilesCode = (String) tmpInnerMap.get(ErtlFunctionalGroupsFinderEvaluationTest.SMILES_CODE_KEY);
            String tmpPseudoSmilesCode = (String) tmpInnerMap.get(ErtlFunctionalGroupsFinderEvaluationTest.PSEUDO_SMILES_CODE_KEY);
            //Writing the record for this functional group
            String tmpRecord = tmpHashCode 
                    + ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR 
                    + tmpPseudoSmilesCode 
                    + ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR 
                    + tmpSmilesCode;
            for (String tmpSettingsKey : this.settingsKeysList) {
                if (tmpInnerMap.get(tmpSettingsKey) == null) {
                    tmpInnerMap.put(tmpSettingsKey, 0);
                }
                tmpRecord += ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR 
                        + tmpInnerMap.get(tmpSettingsKey);
            }
            IAtomContainer tmpMoleculeOfOrigin = (IAtomContainer)tmpInnerMap.get(
                    ErtlFunctionalGroupsFinderEvaluationTest.MOLECULE_OF_ORIGIN_KEY);
            String tmpChebiId = tmpMoleculeOfOrigin.getProperty("ChEBI ID");
            String tmpChemblId = tmpMoleculeOfOrigin.getProperty("chembl_id");
            String tmpCdkTitle = tmpMoleculeOfOrigin.getProperty(CDKConstants.TITLE);
            if (tmpChebiId != null) {
                tmpRecord += ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR + tmpChebiId;
            } else if (tmpChemblId != null) {
                tmpRecord += ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR + tmpChemblId;
            } else if (tmpCdkTitle != null) {
                tmpRecord += ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR + tmpCdkTitle;
            } else {
                tmpRecord += ErtlFunctionalGroupsFinderEvaluationTest.OUTPUT_FILE_SEPERATOR 
                        + ErtlFunctionalGroupsFinderEvaluationTest.MOLECULE_OF_ORIGIN_ID_PLACEHOLDER;
            }
            this.dataOutputPrintWriter.println(tmpRecord);
            this.dataOutputPrintWriter.flush();
            tmpFunctionalGroupsIterator.remove();
        }
        this.dataOutputPrintWriter.close();
        this.exceptionsPrintWriter.close();
        this.filteredMoleculesPrintWriter.close();
    }
    
    /**
     * Returns the Pseudo SMILES code of the given molecule. Pseudo atoms are represented by 'R' atoms and aromatic 
     * (C, S, O, P, Se, N) atoms will be marked by *.
     * 
     * @param aMolecule the molecule to be represented by the Pseudo SMILES string
     * @return the Pseudo SMILES representation of the given molecule
     * @throws CDKException if SmilesGenerator.create() throws a CDKException
     * @throws CloneNotSupportedException if IAtomContainer.clone() throws a CloneNotSupportedException when 
     * invoked on aMolecule
     */
    private String getPseudoSmilesCode(IAtomContainer aMolecule) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpMolecule = aMolecule.clone();
        for (IAtom tmpAtom : tmpMolecule.atoms()) {
            if (tmpAtom.isAromatic()) {
                IAtom tmpReplacementAtom = null;
                if (tmpAtom.getSymbol() != null 
                        && this.pseudoSmilesAromaticElementToPlaceholderElementMap.containsKey(tmpAtom.getSymbol())) {
                    tmpReplacementAtom = new Atom(this.pseudoSmilesAromaticElementToPlaceholderElementMap.get(tmpAtom.getSymbol()));
                }
                if (tmpReplacementAtom != null) {
                    Integer tmpImplicitHydrogenCount = tmpAtom.getImplicitHydrogenCount();
                    AtomContainerManipulator.replaceAtomByAtom(tmpMolecule, tmpAtom, tmpReplacementAtom);
                    tmpReplacementAtom.setImplicitHydrogenCount(
                            tmpImplicitHydrogenCount == null ? 0 : tmpImplicitHydrogenCount);
                }
            }
            tmpAtom.setIsAromatic(false);
            if (tmpAtom instanceof IPseudoAtom && "R".equals(((IPseudoAtom)tmpAtom).getLabel())) {  
                //second condition: see creation of R atoms in ErtlFunctionalGroupsFinder
                IAtom tmpReplacementAtom = new Atom(this.pseudoSmilesAromaticElementToPlaceholderElementMap.get("R"));
                Integer tmpImplicitHydrogenCount = tmpAtom.getImplicitHydrogenCount();
                AtomContainerManipulator.replaceAtomByAtom(tmpMolecule, tmpAtom, tmpReplacementAtom);
                tmpReplacementAtom.setImplicitHydrogenCount(
                        tmpImplicitHydrogenCount == null ? 0 : tmpImplicitHydrogenCount);
            }
        }
        for (IBond tmpBond : tmpMolecule.bonds()) {
            tmpBond.setIsAromatic(false);
        }
        String tmpPseudoSmilesCode = this.smilesGenerator.create(tmpMolecule);
        for (String tmpPlaceholderElementSymbol : this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.keySet()) {
            tmpPseudoSmilesCode = tmpPseudoSmilesCode.replaceAll("(\\[" + tmpPlaceholderElementSymbol + "\\])", 
                    this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.get(tmpPlaceholderElementSymbol))
            .replaceAll("(" + tmpPlaceholderElementSymbol + ")", 
                    this.pseudoSmilesPlaceholderElementToPseudoSmilesSymbolMap.get(tmpPlaceholderElementSymbol));
        }
        return tmpPseudoSmilesCode;
    }
    
    /**
     * Logs molecules that are filtered from the SD file to the filtered molecules file with SMILES code, ChEBI name 
     * and id or ChEMBL id or CDK title and why they were filtered.
     * 
     * @param aMolecule the filtered molecule to be logged
     * @param aCounter the number of filtered molecules so far (will be written to file)
     * @param aCause why this molecule was filtered
     */
    private void logFilteredMolecule(IAtomContainer aMolecule, int aCounter, String aCause) {
        if(!this.areFileOperationsActivated) {
            System.out.println("\nFile operations are not activated, invokation of logFilteredMolecule() is therefore not possible.");
            return;
        }
        this.filteredMoleculesPrintWriter.println();
        this.filteredMoleculesPrintWriter.println(aCounter + ". Filtered Molecule");
        try {
            this.filteredMoleculesPrintWriter.println("SMILES code: " + this.smilesGenerator.create(aMolecule));
        } catch (CDKException | NullPointerException anException){
            this.filteredMoleculesPrintWriter.println("SMILES code: " 
                    + ErtlFunctionalGroupsFinderEvaluationTest.SMILES_CODE_PLACEHOLDER);
        }
        String tmpChebiName = aMolecule.getProperty("ChEBI Name");
        if (tmpChebiName != null)
            this.filteredMoleculesPrintWriter.println("ChEBI name: " + tmpChebiName);
        String tmpChebiId = aMolecule.getProperty("ChEBI ID");
        if (tmpChebiId != null)
            this.filteredMoleculesPrintWriter.println("ChEBI ID: " + tmpChebiId);
        String tmpChemblId = aMolecule.getProperty("chembl_id");
        if (tmpChemblId != null)
            this.filteredMoleculesPrintWriter.println("ChEMBL ID: " + tmpChemblId);
        String tmpCdkTitle = aMolecule.getProperty(CDKConstants.TITLE);
        if (tmpCdkTitle != null)
            this.filteredMoleculesPrintWriter.println("CDK title: " + tmpCdkTitle);
        this.filteredMoleculesPrintWriter.println("Cause: " + aCause);
        this.filteredMoleculesPrintWriter.flush();
    }
    
    /**
     * Logs molecules that raised exceptions somewhere in the processing to the exceptions log file with exception 
     * message and stack trace, SMILES code ChEBI name and id or ChEMBL id or CDK title.
     * 
     * @param anException the exception caused by the molecule
     * @param aSettingsKey a string representation of the settings tested in the current iteration, 
     * e.g. the aromaticity model
     * @param aMolecule the exception-causing molecule
     */
    private void logException(Exception anException, String aSettingsKey, IAtomContainer aMolecule) {
        if(!this.areFileOperationsActivated) {
            System.out.println("\nFile operations are not activated, invokation of logException() is therefore not possible.");
            return;
        }
        this.exceptionsCounter++;
        this.exceptionsPrintWriter.println();
        this.exceptionsPrintWriter.println(this.exceptionsCounter + ". " + anException.getClass() + ": " 
                + anException.getLocalizedMessage());
        this.exceptionsPrintWriter.println("Settings key: " + aSettingsKey);
        try {
            this.exceptionsPrintWriter.println("SMILES code: " + this.smilesGenerator.create(aMolecule));
        } catch (CDKException | NullPointerException aNewException){
            this.exceptionsPrintWriter.println("SMILES code: " 
                    + ErtlFunctionalGroupsFinderEvaluationTest.SMILES_CODE_PLACEHOLDER);
        }
        String tmpChebiName = aMolecule.getProperty("ChEBI Name");
        if (tmpChebiName != null)
            this.exceptionsPrintWriter.println("ChEBI name: " + tmpChebiName);
        String tmpChebiId = aMolecule.getProperty("ChEBI ID");
        if (tmpChebiId != null)
            this.exceptionsPrintWriter.println("ChEBI ID: " + tmpChebiId);
        String tmpChemblId = aMolecule.getProperty("chembl_id");
        if (tmpChemblId != null)
            this.exceptionsPrintWriter.println("ChEMBL ID: " + tmpChemblId);
        String tmpCdkTitle = aMolecule.getProperty(CDKConstants.TITLE);
        if (tmpCdkTitle != null)
            this.exceptionsPrintWriter.println("CDK title: " + tmpCdkTitle);
        anException.printStackTrace(this.exceptionsPrintWriter);
        this.exceptionsPrintWriter.flush();
    }
    //</editor-fold>
}

//<editor-fold defaultstate="collapsed" desc="Enum CustomAtomEncoder">
/**
 * Custom Enumeration of atom encoders for seeding atomic hash codes.
 *
 * @author Jonas Schaub
 * @see BasicAtomEncoder
 * @see AtomEncoder
 */
enum CustomAtomEncoder implements AtomEncoder {
    
    /**
     * Encode if an atom is aromatic or not. This specification is necessary to distinguish functional groups with
     * aromatic environments and those without. For example: [H]O[C] and [H]OC* (pseudo SMILES codes) should be
     * assigned different hash codes by the MoleculeHashGenerator.
     *
     * @see IAtom#isAromatic()
     */
    AROMATICITY {
        
        /**
         *{@inheritDoc}
         */
        @Override
        public int encode(IAtom anAtom, IAtomContainer aContainer) {
            return anAtom.isAromatic()? 3 : 2;
        }
    };
}
//</editor-fold>
