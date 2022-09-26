/**
 * Performance test for
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
package de.ibci.ertlfxgroupsfinder.performancetest;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;

/**
 * An application for testing the performance of the ErtlFunctionalGroupsFinder.find() method under parallelization on 
 * multiple threads. 
 * 
 * @author Jonas Schaub
 * @version 1.0.0.1
 */
public class ErtlFunctionalGroupsFinderPerformanceTest {
    
    //<editor-fold defaultstate="collapsed" desc="Private static final constants">
    /**
     * Name of file for logging occurred exceptions
     */
    private static final String EXCEPTIONS_LOG_FILE_NAME = "Exceptions_Log.txt";
    
    /**
     * Name of file for writing results
     */
    private static final String RESULTS_FILE_NAME = "Results.txt";
    
    /**
     * All allowed atomic numbers to pass to the ErtlFunctionalGroupsFinder;
     * String will be split and resulting integers passed to a set
     */
    private static final String NON_METALLIC_ATOMIC_NUMBERS = "1,2,6,7,8,9,10,15,16,17,18,34,35,36,53,54,86";
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * All allowed atomic numbers to pass to the ErtlFunctionalGroupsFinder as a set of integers (will be parsed from
     * NON_METALLIC_ATOMIC_NUMBERS)
     */
    private Set<Integer> nonMetallicAtomicNumbersSet;
    
    /**
     * The working directory (the jar-file's directory)
     */
    private String workingPath;
    
    /**
     * The given number of different threads to use
     */
    private int numberOfThreadsToUse;
    
    /**
     * All molecules loaded from the SD file
     */
    private IAtomContainer[] moleculesArray;
    
    /**
     * The aromaticity model in use
     */
    private Aromaticity aromaticityModel;
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Instantiates and starts the application. It first loads all molecules from a given SD file into memory and then
     * distributes them equally on the given number of different threads to use. It measures the time it takes for all 
     * threads to complete the extraction of functional groups using the ErtlFunctionalGroupsFinder. It exits the system 
     * if an unexpected exception occurs that prevents the application from working, e.g. an IllegalArgumentException 
     * (will be logged to a file, not printed on the console).
     *
     * @param anArgs the command line arguments, anArgs[0] must be the name of the SD file to load (must be located in 
     * the same directory as the application's JAR file) and anArgs[1] must be the number of different threads to use
     * @throws java.io.IOException if the constructor is unable to open a text file for logging occurred exceptions
     */
    public ErtlFunctionalGroupsFinderPerformanceTest(String[] anArgs) throws IOException {
        this.workingPath = (new File("").getAbsoluteFile().getAbsolutePath()) + File.separator;
        LocalDateTime tmpDateTime = LocalDateTime.now();
        String tmpTimeStamp = tmpDateTime.format(DateTimeFormatter.ofPattern("uuuu_MM_dd_HH_mm"));
        File tmpExceptionsLogFile = new File(this.workingPath 
                + ErtlFunctionalGroupsFinderPerformanceTest.EXCEPTIONS_LOG_FILE_NAME);
        FileWriter tmpExceptionsLogFileWriter = new FileWriter(tmpExceptionsLogFile, true);
        PrintWriter tmpExceptionsPrintWriter = new PrintWriter(tmpExceptionsLogFileWriter);
        tmpExceptionsPrintWriter.println("#########################################################################");
        tmpExceptionsPrintWriter.println("Time-stamp: " + tmpTimeStamp);
        tmpExceptionsPrintWriter.println();
        tmpExceptionsPrintWriter.flush();
        tmpExceptionsPrintWriter.close();
        try {
            if (anArgs.length != 2) {
                throw new IllegalArgumentException("Two arguments (a file name and the number of threads to use) are required.");
            }
            this.numberOfThreadsToUse = 0;
            try {
                this.numberOfThreadsToUse = Integer.parseInt(anArgs[1]);
            } catch (NumberFormatException aNumberFormatException) {
                throw new IllegalArgumentException("Argument \"" + anArgs[1] + "\" must be an integer.");
            }
            if (this.numberOfThreadsToUse <= 0) {
                throw new IllegalArgumentException("The number of threads to use must be at least 1.");
            }
            File tmpDBFile = new File(this.workingPath + anArgs[0]);
            FileInputStream tmpDBFileInputStream;
            try {
                tmpDBFileInputStream = new FileInputStream(tmpDBFile);
            } catch (FileNotFoundException | SecurityException anException) {
                throw new IllegalArgumentException("The database file (name) is not valid: " + anException.getMessage());
            }
            CycleFinder tmpCycleFinder = Cycles.or(Cycles.all(), Cycles.vertexShort());
            this.aromaticityModel = new Aromaticity(ElectronDonation.daylight(), tmpCycleFinder);
            String[] tmpMetalNumbersStrings = ErtlFunctionalGroupsFinderPerformanceTest.NON_METALLIC_ATOMIC_NUMBERS.split(",");
            Integer[] tmpMetalNumbersInt = new Integer[tmpMetalNumbersStrings.length];
            for (int i = 0; i < tmpMetalNumbersStrings.length; i++) {
                tmpMetalNumbersInt[i] = Integer.parseInt(tmpMetalNumbersStrings[i]);
            }
            this.nonMetallicAtomicNumbersSet = new HashSet(Arrays.asList(tmpMetalNumbersInt));
            File tmpResultsLogFile = new File(this.workingPath 
                    + ErtlFunctionalGroupsFinderPerformanceTest.RESULTS_FILE_NAME);
            FileWriter tmpResultsLogFileWriter = new FileWriter(tmpResultsLogFile, true);
            PrintWriter tmpResultsPrintWriter = new PrintWriter(tmpResultsLogFileWriter);
            tmpResultsPrintWriter.println("#########################################################################");
            tmpResultsPrintWriter.println("Time-stamp: " + tmpTimeStamp);
            tmpResultsPrintWriter.println();
            tmpResultsPrintWriter.println("Application initialized. Loading database file named " + anArgs[0] + ".");
            tmpResultsPrintWriter.flush();
            System.out.println("\nApplication initialized. Loading database file named " + anArgs[0] + ".");
            IteratingSDFReader tmpDBReader = new IteratingSDFReader(tmpDBFileInputStream,
                    SilentChemObjectBuilder.getInstance(),
                    true);
            List<IAtomContainer> tmpMoleculesList = new LinkedList<>();
            while (tmpDBReader.hasNext()) {
                try {
                    IAtomContainer tmpMolecule = (IAtomContainer) tmpDBReader.next();
                    tmpMolecule = this.applyFiltersAndPreprocessing(tmpMolecule);
                    tmpMoleculesList.add(tmpMolecule);
                } catch (Exception anException) {
                    /*If an IllegalArgumentException is thrown in applyFiltersAndPreprocessing (meaning that the molecule 
                    should be filtered) the molecule is skipped by catching this exception*/
                }
            }
            try {
                tmpDBReader.close();
            } catch (IOException ex) {
                this.appendToLogfile(ex);
            }
            long tmpSeed = System.nanoTime();
            Collections.shuffle(tmpMoleculesList, new Random(tmpSeed));
            this.moleculesArray = new IAtomContainer[tmpMoleculesList.size()];
            this.moleculesArray = tmpMoleculesList.toArray(this.moleculesArray);
            tmpResultsPrintWriter.println("\nDone Loading database. Found and processed " + this.moleculesArray.length + " valid molecules.");
            System.out.println("Done Loading database. Found and processed " + this.moleculesArray.length + " valid molecules.");
            tmpResultsPrintWriter.flush();
            final int tmpNumberOfMolecules = this.moleculesArray.length;
            int tmpNumberOfMoleculesPerThread = tmpNumberOfMolecules/this.numberOfThreadsToUse;
            List<ExtractFunctionalGroupsTask> tmpListOfThreads = new LinkedList<>();
            //the modulo of dividing the number of valid molecules by the number of threads is simply discarded
            int tmpLastEndIndex = tmpNumberOfMolecules - tmpNumberOfMolecules % this.numberOfThreadsToUse - 1;
            for (int i = 0; i <= tmpLastEndIndex; i += tmpNumberOfMoleculesPerThread) {
                IAtomContainer[] tmpMoleculesArrayForThread = Arrays.copyOfRange(this.moleculesArray, i, i + tmpNumberOfMoleculesPerThread);
                tmpListOfThreads.add(new ExtractFunctionalGroupsTask(tmpMoleculesArrayForThread));
            }
            ExecutorService executor = Executors.newFixedThreadPool(this.numberOfThreadsToUse);
            List<Future<Integer>> tmpFuturesList = new LinkedList<>();
            long tmpStartTime = System.currentTimeMillis();
            try {
                tmpFuturesList = executor.invokeAll(tmpListOfThreads);
            } catch (Exception ex) {
                this.appendToLogfile(ex);
            }
            long tmpEndTime = System.currentTimeMillis();
            tmpResultsPrintWriter.println("Divided molecules onto " + this.numberOfThreadsToUse
                    + " threads. Extraction of functional groups took " + (tmpEndTime - tmpStartTime) + " ms.");
            System.out.println("Divided molecules onto " + this.numberOfThreadsToUse
                    + " threads. Extraction of functional groups took " + (tmpEndTime - tmpStartTime) + " ms.");
            int tmpExceptionsCounter = 0;
            for (Future<Integer> tmpFuture : tmpFuturesList) {
                tmpExceptionsCounter += tmpFuture.get();
            }
            tmpResultsPrintWriter.println(tmpExceptionsCounter + " molecules produced an exception.");
            tmpResultsPrintWriter.flush();
            executor.shutdown();
            tmpResultsPrintWriter.println();
            tmpResultsPrintWriter.flush();
            tmpResultsPrintWriter.close();
        } catch (Exception anException) {
            this.appendToLogfile(anException);
            anException.printStackTrace(System.err);
            System.exit(1);
        }
    }
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Performs all preprocessing needed for the ErtlFunctionalGroupsFinder and throws an IllegalArgumentException
     * if the given molecule should not be passed on to the find() method (filtering).
     *
     * @throws IllegalArgumentException if the given molecule should be filtered i.e. does not meet the
     * ErtlFunctionalGroupsFinder's requirements
     * @throws CDKException for several causes connected to different CDK functionalities
     */
    private IAtomContainer applyFiltersAndPreprocessing(IAtomContainer aMolecule) throws IllegalArgumentException, CDKException {
        if (aMolecule.getAtomCount() == 0 || aMolecule.getBondCount() == 0) {
            throw new IllegalArgumentException("Molecule must be filtered!");
        }
        aMolecule.removeProperty(CDKConstants.CTAB_SGROUPS);
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
        if (!ConnectivityChecker.isConnected(aMolecule)) {
            IAtomContainerSet tmpFragmentsSet = ConnectivityChecker.partitionIntoMolecules(aMolecule);
            IAtomContainer tmpBiggestFragment = null;
            for (IAtomContainer tmpFragment : tmpFragmentsSet.atomContainers()) {
                if (tmpBiggestFragment == null || tmpBiggestFragment.getAtomCount() < tmpFragment.getAtomCount()) {
                    tmpBiggestFragment = tmpFragment;
                }
            }
            aMolecule = tmpBiggestFragment;
        }
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (!this.nonMetallicAtomicNumbersSet.contains(tmpAtom.getAtomicNumber())) {
                throw new IllegalArgumentException("Molecule must be filtered!");
            }
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
        this.aromaticityModel.apply(aMolecule);
        return aMolecule;
    }
    
    /**
     * Appends the given exception's stack trace to a log file.
     *
     * @param anException the exception to log
     */
    private void appendToLogfile(Exception anException) {
        if (anException == null) {
            return;
        }
        PrintWriter tmpPrintWriter = null;
        try {
            FileWriter tmpFileWriter = new FileWriter(this.workingPath 
                    + ErtlFunctionalGroupsFinderPerformanceTest.EXCEPTIONS_LOG_FILE_NAME, 
                    true);
            tmpPrintWriter = new PrintWriter(tmpFileWriter);
            StringWriter tmpStringWriter = new StringWriter();
            anException.printStackTrace(new PrintWriter(tmpStringWriter));
            String tmpStackTrace = tmpStringWriter.toString();
            tmpPrintWriter.println(tmpStackTrace);
        } catch (IOException anIOException) {
            anIOException.printStackTrace(System.err);
        } finally {
            if (tmpPrintWriter != null) {
                tmpPrintWriter.println();
                tmpPrintWriter.flush();
                tmpPrintWriter.close();
            }
        }
    }
    //</editor-fold>
}
