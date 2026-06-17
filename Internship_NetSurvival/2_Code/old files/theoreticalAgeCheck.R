} else if (age_option == "Luo") {
  # Option Luo: Deterministic continuous distribution based on Luo et al. 2023
  n <- 10000000
  # 1. Calculate exact patient counts per bracket
  n1 <- round(n * 0.462)
  n2 <- round(n * 0.519)
  n3 <- n - n1 - n2 # Forces total to equal exactly n, resolving the 100.1% sum
  
  # 2. Generate perfectly spaced sequences
  age1 <- seq(from = 50, to = 64.999, length.out = n1)
  age2 <- seq(from = 65, to = 84.999, length.out = n2)
  age3 <- seq(from = 85, to = 94.999, length.out = n3)
  
  # 3. Combine and shuffle the dataframe to prevent strict age-sorting artifacts
  age <- sample(c(age1, age2, age3))
  
  hist(age, breaks=100)
  
  mean(age)
  
} else if (age_option == "LuoTrunc") {
  # Option LuoTrunc: Deterministic discrete distribution based on Luo et al. 2023
  
  n1 <- round(n * 0.462)
  n2 <- round(n * 0.519)
  n3 <- n - n1 - n2 
  
  age1 <- trunc(seq(from = 50, to = 64.999, length.out = n1))
  age2 <- trunc(seq(from = 65, to = 84.999, length.out = n2))
  age3 <- trunc(seq(from = 85, to = 94.999, length.out = n3))
  
  age <- sample(c(age1, age2, age3))
  
  hist(age, breaks=100)
  
  mean(age)
  