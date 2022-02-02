## FAQ:

Can I/we use this tool to generate individual centile/percentile scores in our clinic/centre/study?
The online tool provides the opportunity to generate centile scores for dataset outside of the original model (and for those within) and we recommend a minimum sample size of at least 100 individuals without a diagnosis acquired with the same protocol to estimate your sites/clinics parameterisation. It is not (yet) suited to generate centile/percentile scores for individuals without first establishing that site specific baseline and we do not recommend using smaller datasets to do so. 

Can I get individual level confidence intervals on my (newly generated) centile scores?
Currently the bootstrapping procedure required to generate individual level confidence intervals (such as those shown in the manuscript) are too computationally intensive to be run through our online interface without significantly restricting access. Depending on our availability we can assist with validation and bootstrapping of new data to generate individual confidence intervals offline and should this be of interest we do encourage people to reach out.

The app doesn’t work for me. I can see the top plot but not the centiles?
In most cases this is likely caused by a column name mismatch between what the app expects and the data provided. First please check our template to make sure these match. Second, if your (multi-site) dataset includes either no control subjects at each site or does not contain at least 10 control subjects per sex per site the app will not compute centiles. 

Can I access your data? Is it openly available?
Unfortunately due to various agreements signed as part of this project and the fact that our database includes unpublished proprietary imaging datasets we are not able to share all data openly. The website and paper contain references and links to all datasets included where you could obtain or request access to the listed datasets. Strictly for non-commercial use all models are made available through GitHub alongside the code-base and tutorials.

Can I contribute data?
We are always keen to improve and update the models and increase the representativeness of said models. Do reach out to us directly to discuss the best ways to organise this.

I have run the app and see the centiles, but when I download my centiles there are only centile for a single phenotype?
We currently have the app setup to compute centiles when a certain feature/phenotype is selected to avoid unnecessarily using resources (i.e., not every centile feature might be needed) to compute centile that are not needed. Selecting a different feature/phenotype and then downloading the output again will include only the selected phenotype. 

If I have multi-site data does the app harmonise this automatically? 
The app currently runs site-specific GAMLSS models, however as we note in our accompanying manuscript this harmonisation is optimal mainly for sites that have at least 100 or more subjects per site. Depending on your situation it may be appropropriate to assess site effects first and possibly harmonise data prior to generating centiles on brainchart.io (for example by using an appropriate ComBAT harmonisation strategy first).

Will the models be updated regularly?
We are continuously receiving more data to update our models and given the computational burden to generate models and validate them we plan to update them at least every 6 months.

Does the app support other brain imaging features, such as for example regional features?
We are always working on expanding the set of available features but have no fixed timeline on when these might be available.

I have an idea about how to use this for a new question/project, can I use your data or do you want to collaborate on it?
We are always open to collaborate though we should note that our time and resources are obviously limited. We strongly encourage you to reach out directly about this and to discuss possibilities and/or make connections to other groups in the consortium that may be able to help or collaborate.
