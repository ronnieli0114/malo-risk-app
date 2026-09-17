### Compute the CLivD score from a data frame

compute_clivd <- function(df) {
  age      <- df$age
  whr      <- df$waist_hip_ratio
  alcohol  <- pmin(df$alcohol_grams_week / 10, 50)
  ggt      <- pmin(df$ggt, 200)
  diabetes <- as.numeric(as.character(df$has_t2dm))
  smoking  <- as.numeric(as.character(df$smoking_binary))
  sex      <- df$sex
  
  rcs_alcohol <- function(x) {
    0.19860813   * x +
      -0.0082096868 * pmax(x - 0.1, 0)^3 +
      0.010575035  * pmax(x - 1,   0)^3 +
      -0.002004756  * pmax(x - 3,   0)^3 +
      -0.00033998925* pmax(x - 9,   0)^3 +
      -2.0602882e-5 * pmax(x - 33,  0)^3
  }
  
  -6.7922721 +
    0.044744302  * age +
    0.32961593   * (whr * 10) +
    rcs_alcohol(alcohol) +
    0.011813962  * ggt +
    0.18721469   * (sex == 2) +
    0.55249734   * (diabetes == 1) +
    0.74679941   * (smoking == 1) +
    0.0054325769 * ggt * (sex == 2) +
    -0.64903176   * (sex == 2) * (smoking == 1)
}


compute_clivd_nonlab <- function(df) {
  age      <- df$age
  whr      <- df$waist_hip_ratio
  alcohol  <- pmin(df$alcohol_grams_week / 10, 50)
  diabetes <- as.numeric(as.character(df$has_t2dm))
  smoking  <- as.numeric(as.character(df$smoking_binary))
  sex      <- df$sex
  
  rcs_alcohol <- function(x) {
    0.19222894     * x +
      -0.00015029544 * pmax(x - 0.1, 0)^3 +
      -0.0021265611  * pmax(x - 1,   0)^3 +
      0.0029832769   * pmax(x - 3,   0)^3 +
      -0.00068765143 * pmax(x - 9,   0)^3 +
      -1.8769011e-5  * pmax(x - 33,  0)^3
  }
  
  -8.0940103 +
    0.044177151  * age +
    0.48927753   * (whr * 10) +
    rcs_alcohol(alcohol) +
    0.69669285   * (diabetes == 1) +
    0.75968055   * (smoking == 1) +
    0.63248362   * (sex == 2) +
    -0.59146649   * (sex == 2) * (smoking == 1)
}

