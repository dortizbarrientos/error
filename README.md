# Uncertainty Propagation via Delta Method

**1. Introduction: Why Predict Evolution with Uncertainty?**

Predicting how populations might evolve in response to environmental changes (like climate change) or selection pressures (like in breeding programs) is a major goal in modern biology. Quantitative genetics provides tools for this, such as Lande's equation ($\Delta \mathbf{z} = \mathbf{G}\boldsymbol{\beta}$), which predicts the change in average trait values ($\Delta \mathbf{z}$) based on genetic architecture ($\mathbf{G}$) and selection ($\boldsymbol{\beta}$).

However, biological measurements are never perfectly precise. We always have some uncertainty about our estimates of genetic variance/covariance (in $\mathbf{G}$) and the strength/direction of selection (in $\boldsymbol{\beta}$). It's crucial to understand how this uncertainty in our input 'ingredients' affects the confidence we have in our final evolutionary prediction ($\Delta \mathbf{z}$). If our input estimates are a bit 'wobbly', how much does our prediction wobble?

This R script provides a practical demonstration of one method – the **multivariate delta method** – used to *approximate* how uncertainty propagates from input parameters to the final prediction. It uses the specific numerical worked example presented in Supplementary Material S3 of the main review paper ["Predictive Evolutionary Genomics: A Probabilistic Framework for Forecasting Adaptation"].

**2. The Biological Scenario: Predicting Evolution of Two Traits**

The code simulates a simple scenario where we want to predict the evolutionary change per generation for two correlated traits (e.g., flowering time $z_1$ and plant height $z_2$) in a population.

* $\mathbf{G}$ (Genetic Variance-Covariance Matrix): Represents the available additive genetic variation for flowering time and height (diagonal elements) and the genetic correlation between them (off-diagonal elements). It describes the genetic 'raw material' for evolution.
* $\boldsymbol{\beta}$ (Selection Gradient Vector): Represents the strength and direction of selection acting on flowering time ($\beta_1$) and height ($\beta_2$). It tells us which trait combinations are favoured in the current environment.
* $\Delta \mathbf{z}$ (Response Vector): The predicted change in the average flowering time ($\Delta z_1$) and average height ($\Delta z_2$) in one generation due to selection acting on the available genetic variation.

**3. What the R Code Does, Step-by-Step:**

* **Section 0 (Setup):** Loads the `ellipse` package, which is needed for drawing the uncertainty ellipse on the plot. If you don't have it, you might need to run `install.packages("ellipse")` once in your R console.
* **Section 1 (Define Parameters):** Sets up the 'best guess' (mean estimated values) for the $\mathbf{G}$ matrix and $\boldsymbol{\beta}$ vector based on the numbers in the S3 worked example. It then calculates the single best-guess prediction for $\Delta \mathbf{z}$.
* **Section 2 (Define Parameter Uncertainty):** Defines the *variance* associated with each parameter estimate. This variance quantifies how uncertain we are about each element in $\mathbf{G}$ and $\boldsymbol{\beta}$. For simplicity in this example, it assumes the uncertainties in the different parameters are independent (no covariance between parameter estimates), so the parameter uncertainty matrix (`Sigma_theta`) is diagonal.
* **Section 3 (Calculate Jacobian):** This is the core mathematical step for the delta method. The Jacobian matrix ($\mathbf{J}_{\theta}$) measures the *sensitivity* of the prediction ($\Delta z_1$ and $\Delta z_2$) to tiny changes in each input parameter ($G_{11}, G_{12}, G_{22}, \beta_1, \beta_2$). For example, it calculates how much the predicted change in flowering time ($\Delta z_1$) would shift if our estimate of the genetic variance for flowering time ($G_{11}$) was slightly different.
* **Section 4 (Delta Method Calculation):** Performs the matrix multiplication: $\Sigma_{\text{pred}} = \mathbf{J}_{\theta} \Sigma_{\theta} \mathbf{J}_{\theta}^T$. This combines the input parameter uncertainty (`Sigma_theta`) with the prediction's sensitivity to those parameters (`J_theta`) to calculate the final *approximated uncertainty* in the prediction (`Sigma_pred`). `Sigma_pred` is a 2x2 matrix containing the variances of $\Delta z_1$ and $\Delta z_2$ (how much each prediction wobbles) and their covariance (how the wobbles are related). The script prints this matrix and the resulting standard deviations.
* **Section 5 (Visualization):** Creates a plot:
    * The **red dot** shows the single best-guess prediction (mean $\Delta \mathbf{z}$).
    * The **blue dashed ellipse** shows the approximate 95% confidence region around the prediction, calculated using `Sigma_pred`. This ellipse visually represents the uncertainty – the true evolutionary outcome likely falls within this region, considering the uncertainty in our input parameters. The size, shape, and orientation of the ellipse depend on the calculated variances and covariance in `Sigma_pred`.

**4. How to Use the Code:**

1.  Save the code provided in the file `error_propagation_delta_method.R`.
2.  Make sure you have R installed, and ideally an interface like RStudio.
3.  Install the `ellipse` package if you haven't already: open R/RStudio and type `install.packages("ellipse")` in the console.
4.  Open the script in R/RStudio and run it (e.g., by clicking 'Source' or typing `source("error_propagation_delta_method.R")` in the console).
5.  The script will print the calculated `Sigma_pred` matrix and standard deviations to the console and generate the plot (usually displayed in the 'Plots' pane in RStudio or potentially saved to a file depending on your R setup).

**5. Biological Interpretation and Context:**

This script demonstrates how uncertainty in estimating fundamental biological parameters (genetic architecture $\mathbf{G}$, selection $\boldsymbol{\beta}$) translates into uncertainty in our evolutionary predictions ($\Delta \mathbf{z}$). The delta method provides a useful analytical *approximation* for this uncertainty.

The resulting plot with the confidence ellipse clearly shows that our prediction isn't a single precise point but rather a region of possibilities. Understanding the size and shape of this region is crucial for making informed decisions. For example, in conservation, a wide uncertainty ellipse might mean we are less confident about whether a population will adapt successfully, potentially requiring more cautious management action. In breeding, it might affect our confidence in the predicted genetic gain from a selection strategy.

While the delta method is helpful, remember it's an approximation best suited for situations where uncertainties are relatively small and relationships are roughly linear. Full Bayesian analyses, as discussed in the main review, provide a more comprehensive way to handle complex uncertainties by simulating directly from the probability distributions of parameters, but they are computationally more intensive. This script serves as a pedagogical tool to illustrate the fundamental concept of error propagation in evolutionary forecasting.

