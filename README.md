# Readme

This a simulation study conducted in order to explain why bifactor models fit significantly better than higher order factor models.

The aim is to test the suggestion by Kan et al (2020) that if a network model is the actual data generating mechanism, a bifactor model will likely outperform the higher order factor model in summarizing the data.

To this end we first simulate data according three models:
- A (second order) g model                                                    
- A bifactor model 
- A network model extracted from the US standardization sample and refitted (confirmatory) on the German standardization sample correlation matrix.

Next, we fit all three models on all these data. This in order to show that the AIC and BIC would pick the true model, if the true model would be in the set of models considered.

Next, we show the results in case the network model would not be considered and the competition is between the higher order factor model and bifactor model. In that case the bifactor indeed beats the higher order factor model.

As, mentioned, the data concerned the WAIS-IV German validation sample data (correlation matrix). With thanks to Jens Lange for providing these data.

## References

[Kan, K.J., de Jonge, H., van der Maas, H.L.J., Levine, S.Z., & Epskamp, S. (2020). How to Compare Latent Factor Models and Psychometric Network Models. *Journal of Intelligence*.](https://github.com/KJKan/mcfarland/blob/master/jintelligence-844295.pdf)

