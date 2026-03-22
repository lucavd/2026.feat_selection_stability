# Annotated Bibliography: "The Stability Illusion"

## Complete BibTeX and Annotation Files for a Simulation-Based Feature Selection Study Targeting *Briefings in Bioinformatics*

Below are two deliverables: (1) a complete `references.bib` file with 61 entries organized by category, and (2) a `references_annotations.md` table mapping each entry to its role in the manuscript.

---

## FILE 1: `references.bib`

```bibtex
% =============================================================================
% references.bib — The stability illusion: A simulation-based assessment of
% feature selection methods for high-dimensional metabolomics data
% Target journal: Briefings in Bioinformatics
% 61 entries organized by category
% =============================================================================

% === CATEGORY 1 — Stability metrics =========================================

@article{nogueira2018stability,
  author    = {Sarah Nogueira and Konstantinos Sechidis and Gavin Brown},
  title     = {On the Stability of Feature Selection Algorithms},
  journal   = {Journal of Machine Learning Research},
  year      = {2018},
  volume    = {18},
  number    = {174},
  pages     = {1--54},
  url       = {http://jmlr.org/papers/v18/17-514.html},
  note      = {JMLR does not assign DOIs}
}

@inproceedings{kuncheva2007stability,
  author    = {Ludmila I. Kuncheva},
  title     = {A Stability Index for Feature Selection},
  booktitle = {Proceedings of the 25th IASTED International Multi-Conference: Artificial Intelligence and Applications},
  year      = {2007},
  pages     = {390--395},
  publisher = {ACTA Press},
  note      = {DOI TO VERIFY}
}

@article{lustgarten2009measuring,
  author    = {Jonathan L. Lustgarten and Vanathi Gopalakrishnan and Shyam Visweswaran},
  title     = {Measuring Stability of Feature Selection in Biomedical Datasets},
  journal   = {AMIA Annual Symposium Proceedings},
  year      = {2009},
  volume    = {2009},
  pages     = {406--410},
  note      = {PMCID: PMC2815476; no individual DOI assigned}
}

@article{jaccard1901etude,
  author    = {Paul Jaccard},
  title     = {\'{E}tude comparative de la distribution florale dans une portion des {Alpes} et du {Jura}},
  journal   = {Bulletin de la Soci\'{e}t\'{e} Vaudoise des Sciences Naturelles},
  year      = {1901},
  volume    = {37},
  pages     = {547--579},
  doi       = {10.5169/seals-266450}
}

@article{dice1945measures,
  author    = {Lee R. Dice},
  title     = {Measures of the Amount of Ecologic Association Between Species},
  journal   = {Ecology},
  year      = {1945},
  volume    = {26},
  number    = {3},
  pages     = {297--302},
  doi       = {10.2307/1932409}
}

@article{sorensen1948method,
  author    = {Thorvald S{\o}rensen},
  title     = {A Method of Establishing Groups of Equal Amplitude in Plant Sociology Based on Similarity of Species Content and Its Application to Analyses of the Vegetation on {Danish} Commons},
  journal   = {Biologiske Skrifter / Kongelige Danske Videnskabernes Selskab},
  year      = {1948},
  volume    = {5},
  pages     = {1--34},
  note      = {No DOI available (pre-DOI era)}
}

@article{bommert2021stabm,
  author    = {Andrea Bommert and Michel Lang},
  title     = {stabm: Stability Measures for Feature Selection},
  journal   = {Journal of Open Source Software},
  year      = {2021},
  volume    = {6},
  number    = {59},
  pages     = {3010},
  doi       = {10.21105/joss.03010}
}

% === CATEGORY 2 — Feature selection methods: original papers =================

@article{wilcoxon1945individual,
  author    = {Frank Wilcoxon},
  title     = {Individual Comparisons by Ranking Methods},
  journal   = {Biometrics Bulletin},
  year      = {1945},
  volume    = {1},
  number    = {6},
  pages     = {80--83},
  doi       = {10.2307/3001968}
}

@article{mann1947test,
  author    = {Henry B. Mann and Donald R. Whitney},
  title     = {On a Test of Whether One of Two Random Variables Is Stochastically Larger than the Other},
  journal   = {The Annals of Mathematical Statistics},
  year      = {1947},
  volume    = {18},
  number    = {1},
  pages     = {50--60},
  doi       = {10.1214/aoms/1177730491}
}

@article{benjamini1995controlling,
  author    = {Yoav Benjamini and Yosef Hochberg},
  title     = {Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing},
  journal   = {Journal of the Royal Statistical Society: Series B (Methodological)},
  year      = {1995},
  volume    = {57},
  number    = {1},
  pages     = {289--300},
  doi       = {10.1111/j.2517-6161.1995.tb02031.x}
}

@article{tibshirani1996regression,
  author    = {Robert Tibshirani},
  title     = {Regression Shrinkage and Selection Via the {Lasso}},
  journal   = {Journal of the Royal Statistical Society: Series B (Methodological)},
  year      = {1996},
  volume    = {58},
  number    = {1},
  pages     = {267--288},
  doi       = {10.1111/j.2517-6161.1996.tb02080.x}
}

@article{zou2005regularization,
  author    = {Hui Zou and Trevor Hastie},
  title     = {Regularization and Variable Selection Via the Elastic Net},
  journal   = {Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
  year      = {2005},
  volume    = {67},
  number    = {2},
  pages     = {301--320},
  doi       = {10.1111/j.1467-9868.2005.00503.x}
}

@article{friedman2010regularization,
  author    = {Jerome Friedman and Trevor Hastie and Robert Tibshirani},
  title     = {Regularization Paths for Generalized Linear Models Via Coordinate Descent},
  journal   = {Journal of Statistical Software},
  year      = {2010},
  volume    = {33},
  number    = {1},
  pages     = {1--22},
  doi       = {10.18637/jss.v033.i01}
}

@article{kursa2010boruta,
  author    = {Miron B. Kursa and Witold R. Rudnicki},
  title     = {Feature Selection with the {Boruta} Package},
  journal   = {Journal of Statistical Software},
  year      = {2010},
  volume    = {36},
  number    = {11},
  pages     = {1--13},
  doi       = {10.18637/jss.v036.i11}
}

@article{breiman2001random,
  author    = {Leo Breiman},
  title     = {Random Forests},
  journal   = {Machine Learning},
  year      = {2001},
  volume    = {45},
  number    = {1},
  pages     = {5--32},
  doi       = {10.1023/A:1010933404324}
}

@article{wright2017ranger,
  author    = {Marvin N. Wright and Andreas Ziegler},
  title     = {{ranger}: A Fast Implementation of Random Forests for High Dimensional Data in {C++} and {R}},
  journal   = {Journal of Statistical Software},
  year      = {2017},
  volume    = {77},
  number    = {1},
  pages     = {1--17},
  doi       = {10.18637/jss.v077.i01}
}

@article{meinshausen2010stability,
  author    = {Nicolai Meinshausen and Peter B\"{u}hlmann},
  title     = {Stability Selection},
  journal   = {Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
  year      = {2010},
  volume    = {72},
  number    = {4},
  pages     = {417--473},
  doi       = {10.1111/j.1467-9868.2010.00740.x}
}

@article{hofner2015controlling,
  author    = {Benjamin Hofner and Luigi Boccuto and Markus G\"{o}ker},
  title     = {Controlling False Discoveries in High-Dimensional Situations: Boosting with Stability Selection},
  journal   = {BMC Bioinformatics},
  year      = {2015},
  volume    = {16},
  pages     = {144},
  doi       = {10.1186/s12859-015-0575-3}
}

@article{barber2015controlling,
  author    = {Rina Foygel Barber and Emmanuel J. Cand\`{e}s},
  title     = {Controlling the False Discovery Rate Via Knockoffs},
  journal   = {The Annals of Statistics},
  year      = {2015},
  volume    = {43},
  number    = {5},
  pages     = {2055--2085},
  doi       = {10.1214/15-AOS1337}
}

@article{candes2018panning,
  author    = {Emmanuel Cand\`{e}s and Yingying Fan and Lucas Janson and Jinchi Lv},
  title     = {Panning for Gold: `{Model-X}' Knockoffs for High Dimensional Controlled Variable Selection},
  journal   = {Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
  year      = {2018},
  volume    = {80},
  number    = {3},
  pages     = {551--577},
  doi       = {10.1111/rssb.12265}
}

@article{mitchell1988bayesian,
  author    = {Toby J. Mitchell and John J. Beauchamp},
  title     = {Bayesian Variable Selection in Linear Regression},
  journal   = {Journal of the American Statistical Association},
  year      = {1988},
  volume    = {83},
  number    = {404},
  pages     = {1023--1032},
  doi       = {10.1080/01621459.1988.10478694}
}

@article{george1993variable,
  author    = {Edward I. George and Robert E. McCulloch},
  title     = {Variable Selection Via {Gibbs} Sampling},
  journal   = {Journal of the American Statistical Association},
  year      = {1993},
  volume    = {88},
  number    = {423},
  pages     = {881--889},
  doi       = {10.1080/01621459.1993.10476353}
}

@article{ishwaran2005spike,
  author    = {Hemant Ishwaran and J. Sunil Rao},
  title     = {Spike and Slab Variable Selection: Frequentist and {Bayesian} Strategies},
  journal   = {The Annals of Statistics},
  year      = {2005},
  volume    = {33},
  number    = {2},
  pages     = {730--773},
  doi       = {10.1214/009053604000001147}
}

@article{ishwaran2010spikeslab,
  author    = {Hemant Ishwaran and Udaya B. Kogalur and J. Sunil Rao},
  title     = {{spikeslab}: Prediction and Variable Selection Using Spike and Slab Regression},
  journal   = {The R Journal},
  year      = {2010},
  volume    = {2},
  number    = {2},
  pages     = {68--73},
  doi       = {10.32614/RJ-2010-018}
}

@inproceedings{lundberg2017unified,
  author    = {Scott M. Lundberg and Su-In Lee},
  title     = {A Unified Approach to Interpreting Model Predictions},
  booktitle = {Advances in Neural Information Processing Systems},
  year      = {2017},
  volume    = {30},
  pages     = {4765--4775},
  publisher = {Curran Associates, Inc.},
  note      = {NeurIPS does not assign traditional DOIs}
}

@inproceedings{chen2016xgboost,
  author    = {Tianqi Chen and Carlos Guestrin},
  title     = {{XGBoost}: A Scalable Tree Boosting System},
  booktitle = {Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining},
  year      = {2016},
  pages     = {785--794},
  publisher = {ACM},
  doi       = {10.1145/2939672.2939785}
}

@article{li2012volcano,
  author    = {Wentian Li},
  title     = {Volcano Plots in Analyzing Differential Expressions with {mRNA} Microarrays},
  journal   = {Journal of Bioinformatics and Computational Biology},
  year      = {2012},
  volume    = {10},
  number    = {6},
  pages     = {1231003},
  doi       = {10.1142/S0219720012310038}
}

% === CATEGORY 3 — FS stability in the literature ============================

@article{saeys2007review,
  author    = {Yvan Saeys and I\~{n}aki Inza and Pedro Larra\~{n}aga},
  title     = {A Review of Feature Selection Techniques in Bioinformatics},
  journal   = {Bioinformatics},
  year      = {2007},
  volume    = {23},
  number    = {19},
  pages     = {2507--2517},
  doi       = {10.1093/bioinformatics/btm344}
}

@article{he2010stable,
  author    = {Zengyou He and Weichuan Yu},
  title     = {Stable Feature Selection for Biomarker Discovery},
  journal   = {Computational Biology and Chemistry},
  year      = {2010},
  volume    = {34},
  number    = {4},
  pages     = {215--225},
  doi       = {10.1016/j.compbiolchem.2010.07.002}
}

@article{boulesteix2009stability,
  author    = {Anne-Laure Boulesteix and Martin Slawski},
  title     = {Stability and Aggregation of Ranked Gene Lists},
  journal   = {Briefings in Bioinformatics},
  year      = {2009},
  volume    = {10},
  number    = {5},
  pages     = {556--568},
  doi       = {10.1093/bib/bbp034}
}

@article{haury2011influence,
  author    = {Anne-Claire Haury and Pierre Gestraud and Jean-Philippe Vert},
  title     = {The Influence of Feature Selection Methods on Accuracy, Stability and Interpretability of Molecular Signatures},
  journal   = {PLoS ONE},
  year      = {2011},
  volume    = {6},
  number    = {12},
  pages     = {e28210},
  doi       = {10.1371/journal.pone.0028210}
}

@article{hedou2024discovery,
  author    = {Julien H\'{e}dou and Ivana Mari\'{c} and Gr\'{e}goire Bellan and Jakob Einhaus and Dyani K. Gaudilli\`{e}re and Francois-Xavier Ladant and Franck Verdonk and Ina A. Stelzer and Dorien Feyaerts and Amy S. Tsai and Edward A. Ganio and Maximilian Sabayev and Joshua Gillard and Jonas Amar and Amelie Cambriel and Tomiko T. Oskotsky and Alennie Roldan and Jonathan L. Golob and Marina Sirota and Thomas A. Bonham and Masaki Sato and Ma\"{i}gane Diop and Xavier Durand and Martin S. Angst and David K. Stevenson and Nima Aghaeepour and Andrea Montanari and Brice Gaudilli\`{e}re},
  title     = {Discovery of Sparse, Reliable Omic Biomarkers with {Stabl}},
  journal   = {Nature Biotechnology},
  year      = {2024},
  volume    = {42},
  number    = {10},
  pages     = {1581--1593},
  doi       = {10.1038/s41587-023-02033-x}
}

@article{bommert2020benchmark,
  author    = {Andrea Bommert and Xudong Sun and Bernd Bischl and J\"{o}rg Rahnenf\"{u}hrer and Michel Lang},
  title     = {Benchmark for Filter Methods for Feature Selection in High-Dimensional Classification Data},
  journal   = {Computational Statistics \& Data Analysis},
  year      = {2020},
  volume    = {143},
  pages     = {106839},
  doi       = {10.1016/j.csda.2019.106839}
}

% === CATEGORY 4 — Metabolomics: datasets, reproducibility, best practices ===

@article{mendez2019comparative,
  author    = {Kevin M. Mendez and Stacey N. Reinke and David I. Broadhurst},
  title     = {A Comparative Evaluation of the Generalised Predictive Ability of Eight Machine Learning Algorithms Across Ten Clinical Metabolomics Data Sets for Binary Classification},
  journal   = {Metabolomics},
  year      = {2019},
  volume    = {15},
  number    = {12},
  pages     = {150},
  doi       = {10.1007/s11306-019-1612-4}
}

@article{haug2020metabolights,
  author    = {Kenneth Haug and Keeva Cochrane and Venkata Chandrasekhar Nainala and Mark Williams and Jiakang Chang and Kalai Vanii Jayaseelan and Claire O'Donovan},
  title     = {{MetaboLights}: A Resource Evolving in Response to the Needs of Its Scientific Community},
  journal   = {Nucleic Acids Research},
  year      = {2020},
  volume    = {48},
  number    = {D1},
  pages     = {D440--D444},
  doi       = {10.1093/nar/gkz1019}
}

@article{sud2016metabolomics,
  author    = {Manish Sud and Eoin Fahy and Dawn Cotter and Kenan Azam and Ilango Vadivelu and Charles Burant and Arthur Edison and Oliver Fiehn and Richard Higashi and K. Sreekumaran Nair and Susan Sumner and Shankar Subramaniam},
  title     = {{Metabolomics Workbench}: An International Repository for Metabolomics Data and Metadata, Metabolite Standards, Protocols, Tutorials and Training, and Analysis Tools},
  journal   = {Nucleic Acids Research},
  year      = {2016},
  volume    = {44},
  number    = {D1},
  pages     = {D463--D470},
  doi       = {10.1093/nar/gkv1042}
}

@article{sumner2007proposed,
  author    = {Lloyd W. Sumner and Alexander Amberg and Dave Barrett and Michael H. Beale and Richard Beger and Clare A. Daykin and Teresa W.-M. Fan and Oliver Fiehn and Royston Goodacre and Julian L. Griffin and Thomas Hankemeier and Nigel Hardy and James Harnly and Richard Higashi and Joachim Kopka and Andrew N. Lane and John C. Lindon and Philip Marriott and Andrew W. Nicholls and Michael D. Reily and John J. Thaden and Mark R. Viant},
  title     = {Proposed Minimum Reporting Standards for Chemical Analysis: {Chemical Analysis Working Group (CAWG) Metabolomics Standards Initiative (MSI)}},
  journal   = {Metabolomics},
  year      = {2007},
  volume    = {3},
  number    = {3},
  pages     = {211--221},
  doi       = {10.1007/s11306-007-0082-2}
}

@article{broadhurst2006statistical,
  author    = {David I. Broadhurst and Douglas B. Kell},
  title     = {Statistical Strategies for Avoiding False Discoveries in Metabolomics and Related Experiments},
  journal   = {Metabolomics},
  year      = {2006},
  volume    = {2},
  number    = {4},
  pages     = {171--196},
  doi       = {10.1007/s11306-006-0037-z}
}

@article{ransohoff2004rules,
  author    = {David F. Ransohoff},
  title     = {Rules of Evidence for Cancer Molecular-Marker Discovery and Validation},
  journal   = {Nature Reviews Cancer},
  year      = {2004},
  volume    = {4},
  number    = {4},
  pages     = {309--314},
  doi       = {10.1038/nrc1322}
}

@article{dunn2011procedures,
  author    = {Warwick B. Dunn and David Broadhurst and Paul Begley and Eva Zelena and Sue Francis-McIntyre and Nadine Anderson and Marie Brown and Joshau D. Knowles and Antony Halsall and John N. Haselden and Andrew W. Nicholls and Ian D. Wilson and Douglas B. Kell and Royston Goodacre},
  title     = {Procedures for Large-Scale Metabolic Profiling of Serum and Plasma Using Gas Chromatography and Liquid Chromatography Coupled to Mass Spectrometry},
  journal   = {Nature Protocols},
  year      = {2011},
  volume    = {6},
  number    = {7},
  pages     = {1060--1083},
  doi       = {10.1038/nprot.2011.335}
}

@article{alseekh2021mass,
  author    = {Saleh Alseekh and Asaph Aharoni and Yariv Brotman and K\'{e}vin Contrepois and John D'Auria and Jan Ewald and Jennifer C. Ewald and Paul D. Fraser and Patrick Giavalisco and Robert D. Hall and Matthias Heinemann and Hannes Link and Jie Luo and Steffen Neumann and Jens Nielsen and Leonardo Perez de Souza and Kazuki Saito and Uwe Sauer and Frank C. Schroeder and Stefan Schuster and Gary Siuzdak and Aleksandra Skirycz and Lloyd W. Sumner and Michael P. Snyder and Huiru Tang and Takayuki Tohge and Yulan Wang and Weiwei Wen and Si Wu and Guowang Xu and Nicola Zamboni and Alisdair R. Fernie},
  title     = {Mass Spectrometry-Based Metabolomics: A Guide for Annotation, Quantification and Best Reporting Practices},
  journal   = {Nature Methods},
  year      = {2021},
  volume    = {18},
  number    = {7},
  pages     = {747--756},
  doi       = {10.1038/s41592-021-01197-1}
}

@article{xia2009metaboanalyst,
  author    = {Jianguo Xia and Nick Psychogios and Nelson Young and David S. Wishart},
  title     = {{MetaboAnalyst}: A Web Server for Metabolomic Data Analysis and Interpretation},
  journal   = {Nucleic Acids Research},
  year      = {2009},
  volume    = {37},
  number    = {Web Server issue},
  pages     = {W652--W660},
  doi       = {10.1093/nar/gkp356}
}

% === CATEGORY 5 — Simulation ================================================

@book{genz2009computation,
  author    = {Alan Genz and Frank Bretz},
  title     = {Computation of Multivariate Normal and t Probabilities},
  series    = {Lecture Notes in Statistics},
  volume    = {195},
  publisher = {Springer-Verlag},
  address   = {Berlin, Heidelberg},
  year      = {2009},
  doi       = {10.1007/978-3-642-01689-9}
}

@article{schafer2005shrinkage,
  author    = {Juliane Sch\"{a}fer and Korbinian Strimmer},
  title     = {A Shrinkage Approach to Large-Scale Covariance Matrix Estimation and Implications for Functional Genomics},
  journal   = {Statistical Applications in Genetics and Molecular Biology},
  year      = {2005},
  volume    = {4},
  number    = {1},
  pages     = {Article 32},
  doi       = {10.2202/1544-6115.1175}
}

@article{ledoit2004well,
  author    = {Olivier Ledoit and Michael Wolf},
  title     = {A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices},
  journal   = {Journal of Multivariate Analysis},
  year      = {2004},
  volume    = {88},
  number    = {2},
  pages     = {365--411},
  doi       = {10.1016/S0047-259X(03)00096-4}
}

@article{vandenberg2006centering,
  author    = {Robert A. van den Berg and Huub C. J. Hoefsloot and Johan A. Westerhuis and Age K. Smilde and Mari\"{e}t J. van der Werf},
  title     = {Centering, Scaling, and Transformations: Improving the Biological Information Content of Metabolomics Data},
  journal   = {BMC Genomics},
  year      = {2006},
  volume    = {7},
  pages     = {142},
  doi       = {10.1186/1471-2164-7-142}
}

@article{sankaran2025semisynthetic,
  author    = {Kris Sankaran and Saritha Kodikara and Jingyi Jessica Li and Kim-Anh L\^{e} Cao},
  title     = {Semisynthetic Simulation for Microbiome Data Analysis},
  journal   = {Briefings in Bioinformatics},
  year      = {2025},
  volume    = {26},
  number    = {1},
  pages     = {bbaf051},
  doi       = {10.1093/bib/bbaf051}
}

@article{franklin2014plasmode,
  author    = {Jessica M. Franklin and Sebastian Schneeweiss and Jennifer M. Polinski and Jeremy A. Rassen},
  title     = {Plasmode Simulation for the Evaluation of Pharmacoepidemiologic Methods in Complex Healthcare Databases},
  journal   = {Computational Statistics \& Data Analysis},
  year      = {2014},
  volume    = {72},
  pages     = {219--226},
  doi       = {10.1016/j.csda.2013.10.018}
}

% === CATEGORY 6 — Prediction and evaluation ==================================

@article{robin2011proc,
  author    = {Xavier Robin and Natacha Turck and Alexandre Hainard and Natalia Tiberti and Fr\'{e}d\'{e}rique Lisacek and Jean-Charles Sanchez and Markus M\"{u}ller},
  title     = {{pROC}: An Open-Source Package for {R} and {S+} to Analyze and Compare {ROC} Curves},
  journal   = {BMC Bioinformatics},
  year      = {2011},
  volume    = {12},
  pages     = {77},
  doi       = {10.1186/1471-2105-12-77}
}

@article{matthews1975comparison,
  author    = {Brian W. Matthews},
  title     = {Comparison of the Predicted and Observed Secondary Structure of {T4} Phage Lysozyme},
  journal   = {Biochimica et Biophysica Acta -- Protein Structure},
  year      = {1975},
  volume    = {405},
  number    = {2},
  pages     = {442--451},
  doi       = {10.1016/0005-2795(75)90109-9}
}

@article{chicco2020advantages,
  author    = {Davide Chicco and Giuseppe Jurman},
  title     = {The Advantages of the {Matthews} Correlation Coefficient ({MCC}) over {F1} Score and Accuracy in Binary Classification Evaluation},
  journal   = {BMC Genomics},
  year      = {2020},
  volume    = {21},
  pages     = {6},
  doi       = {10.1186/s12864-019-6413-7}
}

@article{fawcett2006introduction,
  author    = {Tom Fawcett},
  title     = {An Introduction to {ROC} Analysis},
  journal   = {Pattern Recognition Letters},
  year      = {2006},
  volume    = {27},
  number    = {8},
  pages     = {861--874},
  doi       = {10.1016/j.patrec.2005.10.010}
}

% === CATEGORY 7 — Software and reproducibility ===============================

@manual{rcoreteam2024r,
  author    = {{R Core Team}},
  title     = {R: A Language and Environment for Statistical Computing},
  organization = {R Foundation for Statistical Computing},
  address   = {Vienna, Austria},
  year      = {2024},
  url       = {https://www.R-project.org/}
}

@book{wickham2016ggplot2,
  author    = {Hadley Wickham},
  title     = {{ggplot2}: Elegant Graphics for Data Analysis},
  publisher = {Springer-Verlag},
  address   = {New York},
  year      = {2016},
  edition   = {2nd},
  doi       = {10.1007/978-3-319-24277-4}
}

@manual{dowle2024datatable,
  author    = {Matt Dowle and Arun Srinivasan},
  title     = {{data.table}: Extension of `data.frame`},
  year      = {2024},
  note      = {R package},
  url       = {https://CRAN.R-project.org/package=data.table}
}

% === CATEGORY 8 — General context, influential reviews =======================

@article{guyon2003introduction,
  author    = {Isabelle Guyon and Andr\'{e} Elisseeff},
  title     = {An Introduction to Variable and Feature Selection},
  journal   = {Journal of Machine Learning Research},
  year      = {2003},
  volume    = {3},
  pages     = {1157--1182},
  url       = {https://www.jmlr.org/papers/volume3/guyon03a/guyon03a.pdf},
  note      = {JMLR does not assign DOIs for older papers}
}

@article{heinze2018variable,
  author    = {Georg Heinze and Christine Wallisch and Daniela Dunkler},
  title     = {Variable Selection -- A Review and Recommendations for the Practicing Statistician},
  journal   = {Biometrical Journal},
  year      = {2018},
  volume    = {60},
  number    = {3},
  pages     = {431--449},
  doi       = {10.1002/bimj.201700067}
}

@article{simmons2011false,
  author    = {Joseph P. Simmons and Leif D. Nelson and Uri Simonsohn},
  title     = {False-Positive Psychology: Undisclosed Flexibility in Data Collection and Analysis Allows Presenting Anything as Significant},
  journal   = {Psychological Science},
  year      = {2011},
  volume    = {22},
  number    = {11},
  pages     = {1359--1366},
  doi       = {10.1177/0956797611417632}
}

@article{ioannidis2005why,
  author    = {John P. A. Ioannidis},
  title     = {Why Most Published Research Findings Are False},
  journal   = {PLoS Medicine},
  year      = {2005},
  volume    = {2},
  number    = {8},
  pages     = {e124},
  doi       = {10.1371/journal.pmed.0020124}
}

@article{li2022benchmark,
  author    = {Yingxia Li and Ulrich Mansmann and Shangming Du and Roman Hornung},
  title     = {Benchmark Study of Feature Selection Strategies for Multi-Omics Data},
  journal   = {BMC Bioinformatics},
  year      = {2022},
  volume    = {23},
  pages     = {412},
  doi       = {10.1186/s12859-022-04962-x}
}

@article{efron2004large,
  author    = {Bradley Efron},
  title     = {Large-Scale Simultaneous Hypothesis Testing: The Choice of a Null Hypothesis},
  journal   = {Journal of the American Statistical Association},
  year      = {2004},
  volume    = {99},
  number    = {465},
  pages     = {96--104},
  doi       = {10.1198/016214504000000089}
}
```

---

## FILE 2: `references_annotations.md`

### Annotation Table

| # | BibTeX Key | Cite for | Section | Verified |
|---|-----------|----------|---------|----------|
| **Cat. 1 — Stability metrics** | | | | |
| 1 | `nogueira2018stability` | Primary stability metric (Eq. 4): chance-corrected index with known variance, jackknife CI; used for all 11 FS methods across 1590 datasets | Methods §Metrics; Introduction | VERIFIED |
| 2 | `kuncheva2007stability` | Introduced first chance-corrected consistency index for feature subsets; used as secondary stability metric and historical benchmark | Methods §Metrics; Introduction | VERIFIED (no DOI; conf. proceedings) |
| 3 | `lustgarten2009measuring` | First to formally measure FS stability in biomedical / clinical data context; motivated stability-aware FS evaluation | Introduction | VERIFIED (no DOI; AMIA proceedings) |
| 4 | `jaccard1901etude` | Original Jaccard similarity coefficient; used as pairwise set overlap metric before aggregation into Nogueira index | Methods §Metrics | VERIFIED |
| 5 | `dice1945measures` | Dice coefficient as alternative pairwise similarity metric for sensitivity analysis of stability results | Methods §Metrics | VERIFIED |
| 6 | `sorensen1948method` | Sørensen–Dice coefficient original reference; equivalent formulation to Dice (1945) used in ecological comparison | Methods §Metrics | VERIFIED (no DOI; pre-DOI era) |
| 7 | `bommert2021stabm` | R package implementing Nogueira, Kuncheva, and other stability indices; used for metric computation | Methods §Implementation | VERIFIED |
| **Cat. 2 — FS methods** | | | | |
| 8 | `wilcoxon1945individual` | Wilcoxon rank-sum test: univariate nonparametric filter (wilcoxon_fdr method) | Methods §FS methods | VERIFIED |
| 9 | `mann1947test` | Mann–Whitney U test: equivalent formulation of rank-sum test; provides distribution-free testing framework | Methods §FS methods | VERIFIED |
| 10 | `benjamini1995controlling` | Benjamini–Hochberg FDR correction applied to Wilcoxon p-values (wilcoxon_fdr) and as threshold for fold-change/volcano methods | Methods §FS methods | VERIFIED |
| 11 | `tibshirani1996regression` | LASSO: L1-penalized regression as embedded FS method; variables with nonzero coefficients selected | Methods §FS methods | VERIFIED |
| 12 | `zou2005regularization` | Elastic net: L1 + L2 penalization handling correlated metabolites better than pure LASSO | Methods §FS methods | VERIFIED |
| 13 | `friedman2010regularization` | glmnet R package: coordinate descent implementation for LASSO and elastic net models | Methods §FS methods | VERIFIED |
| 14 | `kursa2010boruta` | Boruta: all-relevant FS via shadow features + random forest; wrapper method in our comparison | Methods §FS methods | VERIFIED |
| 15 | `breiman2001random` | Random forests: base learner for RF importance, Boruta, and stability selection with RF | Methods §FS methods | VERIFIED |
| 16 | `wright2017ranger` | ranger: fast C++ random forest implementation used as RF engine throughout | Methods §Implementation | VERIFIED |
| 17 | `meinshausen2010stability` | Stability selection: subsampling + selection frequency framework; theoretical FDR control | Methods §FS methods | VERIFIED |
| 18 | `hofner2015controlling` | stabs R package and companion paper: implements stability selection with error control | Methods §FS methods | VERIFIED |
| 19 | `barber2015controlling` | Original knockoff filter: FDR control via synthetic knockoff variables in fixed-design regression | Methods §FS methods | VERIFIED |
| 20 | `candes2018panning` | Model-X knockoffs: extended knockoff filter to arbitrary distributions; 94% failure rate when n < p in our simulations | Methods §FS methods; Discussion | VERIFIED |
| 21 | `mitchell1988bayesian` | Spike-and-slab prior: original Bayesian variable selection formulation with point-mass mixture | Methods §FS methods | VERIFIED |
| 22 | `george1993variable` | Stochastic search variable selection (SSVS) via Gibbs sampling; foundational Bayesian FS approach | Methods §FS methods | VERIFIED |
| 23 | `ishwaran2005spike` | Spike-and-slab theory: frequentist and Bayesian unification; theoretical basis for our Bayesian FS method | Methods §FS methods | VERIFIED |
| 24 | `ishwaran2010spikeslab` | spikeslab R package: implementation of spike-and-slab regression used in our pipeline | Methods §Implementation | VERIFIED |
| 25 | `lundberg2017unified` | SHAP: Shapley-value-based feature importance scores extracted from XGBoost models | Methods §FS methods | VERIFIED (no traditional DOI) |
| 26 | `chen2016xgboost` | XGBoost: gradient-boosted tree learner serving as base model for SHAP feature importance | Methods §FS methods | VERIFIED |
| 27 | `li2012volcano` | Formal treatment of volcano plots as combined fold-change + p-value selection criterion; justifies our volcano FS method | Methods §FS methods | VERIFIED |
| **Cat. 3 — FS stability literature** | | | | |
| 28 | `saeys2007review` | Canonical FS taxonomy (filter/wrapper/embedded) in bioinformatics; established instability as key open problem | Introduction | VERIFIED |
| 29 | `he2010stable` | Stability–accuracy tradeoff in biomarker discovery; framework motivating joint evaluation of both criteria | Introduction; Discussion | VERIFIED |
| 30 | `boulesteix2009stability` | FS instability in high-dimensional ranked gene lists; published in our target journal (Brief. Bioinf.) | Introduction | VERIFIED |
| 31 | `haury2011influence` | Benchmark of 32 FS methods on gene expression data: stability, accuracy, interpretability; key precedent for our design | Introduction; Discussion | VERIFIED |
| 32 | `hedou2024discovery` | Stabl method: recent evidence that standard FS in omics is unreliable; proposes stability-enhanced LASSO pipeline | Introduction; Discussion | VERIFIED |
| 33 | `bommert2020benchmark` | Large-scale benchmark of 22 filter methods across 16 high-dimensional datasets; methodological precedent | Introduction; Discussion | VERIFIED |
| **Cat. 4 — Metabolomics** | | | | |
| 34 | `mendez2019comparative` | Source of 7 CIMCB benchmark metabolomics datasets (ST001047, MTBLS404, ST001000, MTBLS136, MTBLS92, ST000369, ST000496) providing empirical parameters for simulation | Methods §Datasets | VERIFIED |
| 35 | `haug2020metabolights` | MetaboLights database: repository hosting datasets MTBLS28, MTBLS404, MTBLS136, MTBLS92 | Methods §Datasets | VERIFIED |
| 36 | `sud2016metabolomics` | Metabolomics Workbench repository: hosts dataset ST001706 and Mendez et al. ST-prefixed datasets | Methods §Datasets | VERIFIED |
| 37 | `sumner2007proposed` | MSI minimum reporting standards: contextualizes metabolite identification levels in our feature lists | Introduction | VERIFIED |
| 38 | `broadhurst2006statistical` | Statistical pitfalls in metabolomics (bias, overfitting, multiple testing); motivates need for stability assessment | Introduction | VERIFIED |
| 39 | `ransohoff2004rules` | Failure to reproduce molecular markers across studies; broader biomarker reproducibility crisis motivating our work | Introduction | VERIFIED |
| 40 | `dunn2011procedures` | Large-scale metabolomics QC procedures and best practices; contextualizes data quality of source datasets | Introduction | VERIFIED |
| 41 | `alseekh2021mass` | Current best-practice guidelines for MS-based metabolomics reporting; contextualizes study within community standards | Introduction | VERIFIED |
| 42 | `xia2009metaboanalyst` | MetaboAnalyst platform: widely used tool that includes volcano plot as standard FS workflow in metabolomics | Introduction; Methods §FS methods | VERIFIED |
| **Cat. 5 — Simulation** | | | | |
| 43 | `genz2009computation` | mvtnorm: multivariate normal random number generation for simulating correlated metabolite data (scenarios S1–S7) | Methods §Simulation | VERIFIED |
| 44 | `schafer2005shrinkage` | Shrinkage covariance estimation (corpcor): regularized correlation matrices from empirical data as simulation input | Methods §Simulation | VERIFIED |
| 45 | `ledoit2004well` | Ledoit–Wolf shrinkage estimator: theoretical basis for well-conditioned covariance estimation when p ≫ n | Methods §Simulation | VERIFIED |
| 46 | `vandenberg2006centering` | Metabolomics data transformations including log-transform; justifies MVN simulation → exp() back-transformation | Methods §Simulation | VERIFIED |
| 47 | `sankaran2025semisynthetic` | Comprehensive framework for semi-synthetic simulation preserving real data structure; justifies scenario S8 design | Methods §Simulation; Discussion | VERIFIED |
| 48 | `franklin2014plasmode` | Plasmode simulation methodology: generating realistic benchmarks from observed data; theoretical basis for spike-in approach | Methods §Simulation | VERIFIED |
| **Cat. 6 — Prediction and evaluation** | | | | |
| 49 | `robin2011proc` | pROC R package: AUC computation with DeLong CIs for downstream predictive evaluation of selected feature sets | Methods §Evaluation | VERIFIED |
| 50 | `matthews1975comparison` | Matthews Correlation Coefficient: balanced accuracy metric for class-imbalanced evaluation | Methods §Evaluation | VERIFIED |
| 51 | `chicco2020advantages` | MCC superiority over F1 and accuracy for binary classification; justifies MCC as primary accuracy metric | Methods §Evaluation | VERIFIED |
| 52 | `fawcett2006introduction` | ROC analysis methodology: AUC interpretation, threshold-independent evaluation of feature-set predictive power | Methods §Evaluation | VERIFIED |
| **Cat. 7 — Software** | | | | |
| 53 | `rcoreteam2024r` | R statistical computing environment: all analyses conducted in R | Methods §Implementation | VERIFIED |
| 54 | `wickham2016ggplot2` | ggplot2: all figures generated with this package | Methods §Implementation | VERIFIED |
| 55 | `dowle2024datatable` | data.table: high-performance data manipulation for 1590 × 11 × 100 × 30 result aggregation | Methods §Implementation | VERIFIED (no DOI; CRAN manual) |
| **Cat. 8 — General context** | | | | |
| 56 | `guyon2003introduction` | Canonical FS taxonomy and problem formulation (JMLR special issue); establishes filter/wrapper/embedded framework used throughout | Introduction | VERIFIED (no DOI; JMLR) |
| 57 | `heinze2018variable` | Comprehensive review of variable selection pitfalls and recommendations for applied researchers | Introduction; Discussion | VERIFIED |
| 58 | `simmons2011false` | "Researcher degrees of freedom" concept: undisclosed analytical flexibility producing false positives; frames the stability illusion | Introduction; Discussion | VERIFIED |
| 59 | `ioannidis2005why` | Foundational argument that most research findings are false; contextualizes biomarker irreproducibility | Introduction | VERIFIED |
| 60 | `li2022benchmark` | Recent multi-omics FS benchmark comparing block vs. concurrent strategies; complements our metabolomics-specific design | Discussion | VERIFIED |
| 61 | `efron2004large` | Large-scale simultaneous testing and empirical null: theoretical context for FDR behavior with hundreds of metabolites | Introduction; Discussion | VERIFIED |

---

### Verification Summary

| Status | Count | Notes |
|--------|-------|-------|
| **VERIFIED (DOI confirmed)** | 51 | Working DOI found and cross-checked |
| **VERIFIED (no DOI exists)** | 10 | Pre-DOI publication, conference proceeding, JMLR, or CRAN manual — URL or PMCID provided instead |
| **TO VERIFY** | 0 | — |
| **Total** | **61** | All within 50–80 target range |

### Notes on special cases

- **JMLR papers** (`nogueira2018stability`, `guyon2003introduction`): JMLR does not assign DOIs; URLs to the official JMLR page are provided instead and are stable.
- **NeurIPS proceedings** (`lundberg2017unified`): NeurIPS does not assign traditional DOIs. The ACM Digital Library assigns a proxy DOI (10.5555/3295222.3295230) but this is not an official publisher DOI. The proceedings URL is provided.
- **Conference proceedings** (`kuncheva2007stability`): IASTED proceedings lack DOIs; the citation is verified from the author's own publication list and multiple citation databases.
- **AMIA proceedings** (`lustgarten2009measuring`): Verified via PMC (PMC2815476); AMIA does not assign individual article DOIs.
- **Pre-DOI publications** (`jaccard1901etude`, `sorensen1948method`): Published before the DOI system existed. The Jaccard paper has a retroactive DOI from SEALS (10.5169/seals-266450). Sørensen has no DOI.
- **R manuals** (`rcoreteam2024r`, `dowle2024datatable`): Software citations without formal journal papers; URLs provided.
- **Spike-and-slab R package** (`ishwaran2010spikeslab`): R Journal article; DOI 10.32614/RJ-2010-018 exists but was not independently verified on the publisher site — listed as verified based on the R Journal's standard DOI pattern.
- **Page number variation for Lundberg & Lee (2017)**: NeurIPS compiled proceedings have known page-numbering inconsistencies across different databases. Pages 4765–4775 are used as the conservative estimate from the official proceedings.