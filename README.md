# A mathematical proof comparing the expectation and variance between two common approaches for parameterizing sex-composition data in fishery stock assessments
## Matthew LH. Cheng, Peter-John F. Hulson, Daniel R. Goethel, Curry J. Cunningham
### Contact: lhcheng@alaska.edu
Two primary methods for parameterizing sex-specific age and length composition data in fishery stock assessments exist, which we refer to as the ‘Across’ and ‘Within’ approaches. When using the ‘Across’ approach, sex-composition data are assumed to arise from a single statistical model that describes the probability of sampling across all ages and sexes in a given year. By contrast, the ‘Within’ approach assumes that sex-composition data arises from several statistical models: sex-specific models that describe the probability of sampling ages within each sex, and an additional model that describes the sex-ratio information of composition data. In this mathematical proof, we derive the statistical properties of both approaches under Multinomial and Dirichlet-Multinomial sampling and show that they produce equivalent model expectations. However, we illustrate that the ‘Within’ approach leads to smaller assumed variances when sampling follows a Dirichlet-Multinomial distribution, because overdispersion acts independently within each sex rather than jointly across sexes. Considering that both approaches yield equivalent model expectations, we recommend the ‘Across’ approach for parameterizing sex-composition data due to its simplicity in implementation, its alignment with most fisheries sampling designs, and its ability to jointly account for overdispersion and sampling correlations across sexes.

