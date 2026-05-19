# **Modify the provided code to generate the data by:**

1. **Using an exponential distribution for the baseline excess hazard**
   - Set **lambda = 0.03**,
   - This value may be varied later.
2. **Age distribution: TWO OPTIONS to test** a) age \~ Uniform\[85, 90\]  
    b) age \~ Normal(65, 10) truncated to \[50, 85\]  
    The **proportion of cancer‑related deaths** obtained under both options will need to be compared (see point 7).
3. **Sex and X (binary variable):**
   - Proportion of women in the dataset: **pFem <- 0.5**
   - Proportion of X (=0) in the dataset: **pX <- 0.5**
4. **Effects of covariates on excess mortality**
   - betaSex  <- –1.5
   - betaX    <- 0.01
   - betaAge  <- 0.02
5. **Apply administrative censoring at 5 years**

---

Once the data have been generated for **one dataset (N = 1) of 10,000 patients (n = 10,000):**

1. **Verification steps**
   - Check the distributions of **age**, **sex**, and **X**,
   - Plot the **estimated net survival** together with the **theoretical net survival** (on the same graph), where the theoretical curve is obtained using the **Kaplan–Meier estimator** applied to the hypothetical‑world data,
   - Explore the results of **relsurv::nessie** on these data  
      (Reminder: an example of how to use this function is provided in the Pavlic & Pohar‑Perme paper in *JSS*.)
2. **Additional calculations**
   - compute the **proportion of deaths due to cancer**,
   - compute the **proportion of deaths due to other causes**.

---