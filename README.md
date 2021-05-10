# GraphSolMol
The Graph Solubility Modelling (GraphSolMol) Repository contains the python code for an undergraduate research project at the University of Nottingham on aqueous solubility modelling using vectorial and graph kernels for SVM. For any queries, please contact either myself (Phil Wroe) or Joe Redshaw, the contact details for both can be found below:

  pcypw2@nottingham.ac.uk
  
  psxjr4@nottingham.ac.uk

The Repository is structured in the following way:

# Datasets
The AqSolDB(1) was used as the primary dataset. Further preprocessing was undertaken to remove extreme datapoints that were not relevant to drug discovery. Compounds with a molecular weight > 1000 were removed, and only compounds within a -11 to +11 logP range were considered. 

A validation dataset was curated from a model by Fioressi et al.(2) was curated using the same algorimth used on the AqSolDB(3) to remove dupliates, validate the solubiltiy datapoints and gorup them based on their standard deviation.  

# Descriptors
The physiochemical descriptors, molecular weight (MolWT), lipophilcity (logP), aromatic proportion (AP) and rotatable bonds (RB) from the ESOL model (Delaney at el.) were used as a benchmark model. 
Secondly, the structual binary descriptors MACCS key fingerprints and Extended Connectivity Fingerprints (ECFP) (r=3) were used. Both the structural and physiochemcal descriptors were built from RDkit. These were calculated an exported to a .csv file prior to model building. I've included the code to calculate these descriptors. 
Finally, The graphical representation was built using Grakel(5). 

# Methods 
Recreating the ESOL MLR model formed part of the preliminary analysis of this project. The aims of this work was to firstly build a benchmark model to compare the SVR results against, but to also aid my decision making in creating a dataset that captured the extremities of chemical space in the context of aqueous solubility.  For more information on the exact methodology, please refer to results tab. 

The primary method used was support vector machines (SVM). The vectorial kernels explored were linear, RBF, polynomial (degree=1,2,3), sigmoid. The computation was done in house using scikit-learn, so applying them is relatively straight forward. The graphical kernels were slightly more technical; The kernel methods assessed here were random walk, shortest path, ordered dag decomposition, neighbourhood hashed, propagation, pyramid match, graph sampling, vertex histogram, Weisfeiler-Lehman subtree and Weisfeiler-Lehman optimal assignment. Computation of kernel matrix was done externally using grakel before being inputted into SVM moddel. 

The SVM hypereparameters were optimised using a grid search 5-fold cross validation method and the model was built on a 80:20 split of the AqSolDB. Coefficient of determination and root mean squared error were calculated as part of the model evaluation. 


# Results
I have attached a word file of my results for the project. As part of the project, I was aiming to have gathered results into indefinite graph kernels, but due to time constraints this had to be suspended. Graphical methods appear to be more accurate at modelling solubility based of the predicted vs experimental plots I have attached which is great! There is still a lot of general noise using graphs, which can be attribued by the nature of the solubiltiy data being very sensitive to experimental conditions. If anyone would like to contribute to this project, either by applying a more accurate dataset or by exploring graph kernels in further detail, please feel free! I have included all my code and data for complete transparency and reproducability so it should be straight forward to gather results. Likewise, I don't doubt the code isn't fully optimised (please bare in mind I'm an undergraduate with limited coding experience), so if you have any issues or have any suggestions on how to improve it, then  please email myself (Phil Wroe) or my PhD mentor (Joe Redshaw). 

To reference this any of the results included in this github, please use the following:

Wroe P., Redshaw J., Hirst J. Evaluation of Vectorial and Graph Kernel Method Performance on Aqueous Solubility Predictions. Nottingham J. Chem. (2021). 

# References
(1) Sorkun, M.C., Khetan, A. & Er, S. AqSolDB, a curated reference set of aqueous solubility and 2D descriptors for a diverse set of compounds. Sci Data 6, 143 ##(2019).

(2)Fioressi S.E., Bacelo D.E., Aranda J.F. & Duchowicz P.R. Prediction of the aqueous solubility of diverse compounds by 2D-QSPR. J. Mol. Liq., 302, 112572 (2020).

(3) https://codeocean.com/capsule/8848590/tree/v1

(4) Delaney J.S. ESOL: Estimating Aqueous Solubility Directly from Molecular Structure. J. Chem. Inf. Model, 44(3), 1000-1005 (2004).

(5)https://ysig.github.io/GraKeL/0.1a8/index.html
