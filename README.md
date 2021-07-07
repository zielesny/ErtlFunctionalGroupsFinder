# ErtlFunctionalGroupsFinder

The algorithm for automated functional groups detection and extraction of organic molecules developed by Peter Ertl is implemented on the basis of the Chemistry Development Kit (CDK).<br/>
(Project team: Sebastian Fritsch, Stefan Neumann, Jonas Schaub, Christoph Steinbeck, Achim Zielesny).

Folder **Basic** contains the basic ErtlFunctionalGroupsFinder code and test code for integration in Java projects (see instructions in *Basic usage instructions.txt*).

Folder **CDK** contains CDK library jar file *cdk-2.2.jar* that ErtlFunctionalGroupsFinder works with.

Folder **Evaluation** contains sample code for evaluation of functional groups with ErtlFunctionalGroupsFinder (see instructions in *Evaluation usage instructions.txt*).

Folder **JUnit 4** contains library jar files for unit testing.

Folder **Performance** contains a jar library for performance tests (see instructions in *Performance usage instructions.txt*).

Folder **Utility** contains code for filtering, preprocessing, and postprocessing routines useful for evaluation of functional groups with ErtlFunctionalGroupsFinder (see instructions in *Utility usage instructions.txt*).

ErtlFunctionalGroupsFinder is described in the scientific literature (see reference below).

# Acknowledgments
The authors like to thank Peter Ertl for describing his algorithm in a way that allowed easy re-implementation. This is not always the case. We also thank him for valuable discussions. We appreciate help from Egon Willighagen and John Mayfield with the CDK integration and from Felix Bänsch for unbiased release testing.

# References
**ErtlFunctionalGroupsFinder**<br/>
* [Fritsch S, Neumann S, Schaub J, Steinbeck C, Zielesny A. ErtlFunctionalGroupsFinder: automated rule-based functional group detection with the Chemistry Development Kit (CDK). J Cheminform. 2019; 11:37](https://doi.org/10.1186/s13321-019-0361-8)

**Ertl algorithm**<br/>
* [Ertl P. An algorithm to identify functional groups in organic molecules. J Cheminform. 2017; 9:36.](https://doi.org/10.1186/s13321-017-0225-z)

**Chemistry Development Kit (CDK)**<br/>
* [Chemistry Development Kit on GitHub](https://cdk.github.io/)<br/>
* [Steinbeck C, Han Y, Kuhn S, Horlacher O, Luttmann E, Willighagen EL. The Chemistry Development Kit (CDK): An Open-Source Java Library for Chemo- and Bioinformatics. J Chem Inform Comput Sci. 2003;43(2):493-500.](https://dx.doi.org/10.1021%2Fci025584y)<br/>
* [Steinbeck C, Hoppe C, Kuhn S, Floris M, Guha R, Willighagen EL. Recent Developments of the Chemistry Development Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics. Curr Pharm Des. 2006; 12(17):2111-2120.](https://doi.org/10.2174/138161206777585274)<br/>
* [May JW and Steinbeck C. Efficient ring perception for the Chemistry Development Kit. J. Cheminform. 2014; 6:3.](https://dx.doi.org/10.1186%2F1758-2946-6-3)<br/>
* [Willighagen EL, Mayfield JW, Alvarsson J, Berg A, Carlsson L, Jeliazkova N, Kuhn S, Pluska T, Rojas-Chertó M, Spjuth O, Torrance G, Evelo CT, Guha R, Steinbeck C, The Chemistry Development Kit (CDK) v2.0: atom typing, depiction, molecular formulas, and substructure searching. J Cheminform. 2017; 9:33.](https://doi.org/10.1186/s13321-017-0220-4)
