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

import java.util.concurrent.Callable;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ErtlFunctionalGroupsFinder;

/**
 * A Callable thread that extracts functional groups from all molecules in a given array using the 
 * ErtlFunctionalGroupsFinder class.
 * 
 * @author Jonas Schaub
 */
public class ExtractFunctionalGroupsTask implements Callable<Integer> {

    private final IAtomContainer[] moleculesArray;
    
    private final ErtlFunctionalGroupsFinder ertlFinder;
    
    /**
     * Instantiates the thread.
     * 
     * @param aListOfMolecules atom containers should meet the ErtlFunctionalGroupsFinder's input specifications but 
     * any occurring exception will be caught
     */
    public ExtractFunctionalGroupsTask(IAtomContainer[] aListOfMolecules) {
        this.moleculesArray = aListOfMolecules;
        this.ertlFinder = new ErtlFunctionalGroupsFinder();
    }
    
    /**
     * Applies the ErtlFunctionalGroupsFinder.find(IAtomContainer container, boolean clone) method on all given 
     * molecules (parameter clone = false) and counts the occurring exceptions.
     * 
     * @return the number of occurred exceptions
     * @throws Exception if unable to compute a result (copied from doc in Callable interface)
     */
    @Override
    public Integer call() throws Exception {
        int tmpExceptionsCounter = 0;
        for (IAtomContainer tmpMolecule : this.moleculesArray) {
            try {
                this.ertlFinder.find(tmpMolecule, false);
            } catch (Exception anException) {
                tmpExceptionsCounter++;
            }
        }
        return tmpExceptionsCounter;
    }
    
}
