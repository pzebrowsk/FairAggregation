##Package FairAggregation contains fair aggregating functions

#' Negative of Gini index
#'
#' @description
#' \code{negative_Gini_index} returns the value of Gini index times -1. The higher value of negative Gini index the lower inequality.
#'
#' @param x A numeric vector.
#'
#' @details
#' \code{x} is expected to be a utility profile, i.e., none of its elements is negative. If \code{x} contains negative values, a warning message is displayed, but the value of the function is returned.
#'
#' Formula for negative Gini index:
#' \deqn{
#' -G(x) = -\frac{\sum_{i=1}^N \sum_{i=1}^N |x_i - x_j|}{2N\sum_{j=1}^N x_j}.
#' }
#'
#' For a proper utility profile \eqn{x}, Gini index \eqn{G} takes value between 0 and \eqn{(1-1/N)}, where \eqn{N} is the length of \eqn{x}. Value 0 corresponds to the state of perfect equality (i.e., all elements of \eqn{x} are equal), while value \eqn{(1-1/N)} is attained in the state of extreme inequality (i.e., all but one element of \eqn{x} are equal to zero).
#'
#' @returns
#' Value of negative Gini index \eqn{-G(x)}.
#' @export
#' negative_Gini_index
#'
#' @examples
#' negative_Gini_index(c(1,2,3,4))
#' negative_Gini_index(c(1,0,0,0))
#' negative_Gini_index(c(1,1,1,1))

negative_Gini_index <- function(x) {
  if (any(x < 0)) {
    warning("Gini coefficient has no interpretation for vectors containing negative values!")
  }

  y <- 0        # variable storing sum of difference within the vector x
  for (u in x) {
    for (v in x) {
      y <- y + abs(u-v)
    }
  }
  return(-(y)/(2.0*length(x)*sum(x)))
}


####............................................................................

#' Negative of generalized entropy index
#'
#' @description
#' \code{negative_generalized_entropy_index} returns the value of generalized index times -1. The higher value of negative generalized entropy index the lower inequality.
#'
#' @param x A numeric vector
#' @param alpha A real number
#'
#' @details
#' \code{x} is expected to be a utility profile, i.e., none of its elements is negative. If \code{x} contains negative values, a warning message is displayed, but the value of the function is returned.
#'
#' Formula for negative generalized entropy index:
#' \deqn{
#' -g_\alpha(x) =
#'  \left\{
#'    \begin{array}{ll}
#'      \frac{1}{N} \sum_{i=1}^N \ln{\frac{x_i}{\mu(x)}} &  \qquad \mathrm{for} \quad \alpha = 0\\
#'      -\frac{1}{N} \sum_{i=1}^N \frac{x_i}{\mu(x)} \ln{\frac{x_i}{\mu(x)}} &  \qquad \mathrm{for} \quad \alpha = 1\\
#'      \frac{-1}{N \alpha (\alpha - 1)} \sum_{i=1}^N \left( (\frac{x_i}{\mu(x)})^\alpha -1 \right) &  \qquad \mathrm{for} \quad \alpha \neq 0, 1
#'    \end{array}
#'  \right.
#' }
#' where \eqn{\mu(x)} is the mean of \code{x} and \eqn{N} is the length of \eqn{x}.
#'
#' The parameter \eqn{\alpha} determines the sensitivity to changes in different segments of the utility profile \eqn{x}: lower values imply higher sensitivity to changes in the lower tail of the utility profile, while higher values imply higher sensitivity to changes in the upper tail.
#'
#' @returns
#' Value of negative generalized entropy index \eqn{-g_\alpha(x)}.
#' @export
#' negative_generalized_entropy_index
#'
#' @examples
#' negative_generalized_entropy_index(c(1,2,3,4))
#' negative_generalized_entropy_index(c(1,2,3,4), alpha = 2)
#'
#'
negative_generalized_entropy_index <- function(x, alpha = 1) {
  if (any(x < 0)) {
    warning("Generalized entropy index has no proper interpretation for vectors containing negative values!")
  }

  res <- 0
  m <- mean(x)
  N <- length(x)

  if (alpha == 0) {
    res <- sum(log(x/m))/N
  }
  else if (alpha == 1) {
    res <- -sum((x/m)*log(x/m))/N
  }
  else {
    res <- -sum((x/m)^alpha - 1)/(N*alpha*(alpha - 1))
  }

  return(res)
}

####............................................................................



#' Negative Atkinson index
#'
#' @description
#' \code{negative_Atkinson_index} returns the value of the Atkinson inequality index times -1. The higher value of negative Atkinson index the lower inequality.
#'
#' @param x A numeric vector
#' @param gamma A non-negative real number
#'
#' @details
#' \code{x} is expected to be a utility profile, i.e., none of its elements is negative. If \code{x} contains negative values, a warning message is displayed, but the value of the function is returned.
#'
#' \code{gamma} must be a non-negative real number. If negative value of \code{gamma} is passed, a warning message is displayed and execution of the function stopped.
#'
#' Formula for negative Atkinson index:
#' \deqn{
#' -I_\gamma(x) =
#'  \left\{
#'    \begin{array}{ll}
#'      \left(\frac{N}{\sum_{i=1}^N x_i}\right) \left( \prod_{i=1}^N x_i \right)^{\frac{1}{N}} - 1 &  \qquad \mathrm{for} \quad \gamma = 1\\
#'      \left(\frac{N}{\sum_{i=1}^N x_i}\right) \left( \frac{\sum_{i=1}^N x_i^{1-\gamma}}{N} \right)^{\frac{1}{1-\gamma}} - 1 &  \qquad \mathrm{for} \quad \gamma \geq 0,  \gamma \neq 1
#'    \end{array}
#'  \right.
#' }
#' where \eqn{N} is the length of \eqn{x}.
#'
#' The parameter \eqn{\gamma} determines the social aversion to inequality: \eqn{\gamma = 0} corresponds to no aversion, \eqn{\gamma = \infty} corresponds to full aversion (no tolerance) to inequality.
#'
#' @returns
#' Value of negative Atkinson index \eqn{-I_\gamma(x)}.
#' @export
#' negative_Atkinson_index
#'
#' @examples
#' negative_Atkinson_index(c(1,2,3,4))
#' negative_Atkinson_index(c(1,2,3,4), alpha = 2)
negative_Atkinson_index <- function(x, gamma = 1) {
  if (any(x < 0)) {
    warning("Atkinson index has no proper interpretation for vectors containing negative values!")
  }

  if (gamma < 0) {
    stop("Atkinson index not defined for gamma < 0.")
  }

  res <- 0
  N <- length(x)

  if (gamma == 1) {
    res <- (N/sum(x))*(prod(x))^(1/N) - 1
  }
  else {
    res <- (N/sum(x))*(sum(x^(1-gamma)/N))^(1/(1-gamma)) - 1
  }

  return(res)
}

####............................................................................


#' Generalized mean
#'
#' @description
#' \code{generalized_mean} computes the generalized mean of a vector.
#'
#' @param x A numeric vector
#' @param p A real number
#'
#' @details
#' \code{x} is expected to be a utility profile, i.e., none of its elements is negative. If \code{x} contains negative values, a warning message is displayed, but the value of the function is returned.
#'
#' Formula for generalized mean:
#' \deqn{
#' g_p(x) =
#'  \left\{
#'    \begin{array}{ll}
#'      \left( \prod_{i=1}^N x_i \right)^{\frac{1}{N}} &  \qquad \mathrm{for} \quad p = 0\\
#'      \left(\frac{1}{N} \sum_{i=1}^N x_i^p \right)^{\frac{1}{p}} &  \qquad \mathrm{for} \quad p \neq 0
#'    \end{array}
#'  \right.
#' }
#' where \eqn{N} is the length of \eqn{x}.
#'
#' The parameter \eqn{p} determines the trade-off between equality (egalitarianism for \eqn{p = -\infty}) and efficiency (utilitarianism for \eqn{p = 1}). For values of \eqn{p > 1} generalized mean is a convex function, attaining the higher values the more unequal is the level of inequality in a utility profile \eqn{x}, provided the total utility \eqn{\sum_{i=1}^N x_i} is kept constant.
#'
#' @return
#' Value of generalized mean \eqn{g_p(x)}.
#' @export
#' generalized_mean
#'
#' @examples
#' generalized_mean(c(1,2,3,4))
#' generalized_mean(c(1,2,3,4), alpha = 1)
#'
#'
generalized_mean <- function(x, p = 0) {
  res <- 0
  N <- length(x)

  if (p == 0) {
    res <- (prod(x))^(1/N)
  }
  else {
    res <- (sum(x.^p)/N)^(1/p)
  }

  return(res)
}

####............................................................................

#
# # Fairness ration requires the whole set U of attainable utility profiles
# fairness_ratio <- function(x, U) {
#
# }

####............................................................................


#' Underachievement function
#'
#' @description
#' \code{underachievement_fun} computes the value of underachievement function of a vector.
#'
#' @param x A numeric vector
#' @param alpha A non-negative real number
#' @param ineq_measure A function (standard deviation by default)
#'
#' @details
#' Formula for underachievement function:
#' \deqn{
#' g_\alpha(x) = \mu(x) - \alpha\rho(x)
#' }
#' where \eqn{\mu(x)} is the mean of \eqn{x} and \eqn{\rho(\cdot)} is an inequality measure.
#'
#' The value of underachievement function is the average utility \eqn{\mu(x)} with penalty for inequality in utility distribution, as measured by the inequality measure \eqn{\rho(x)}. The parameter \eqn{\alpha} is the penalty weight. For no inequality penalty, i.e,  \eqn{\alpha = 0}, the underachievement function is equivalent to utilitarian total utility function. The higher the value of \eqn{\alpha} the larger emphasis on equality.
#'
#' @return
#' Value of underachievement function \eqn{g_\alpha(x)}.
#' @export
#' underachievement_fun
#'
#' @examples
#' underachievement_fun(c(1,2,3,4))
#' underachievement_fun(c(1,2,3,4), alpha = 0.2, ineq_measure = var)
#'
#'
underachievement_fun <- function(x, alpha = 0, ineq_measure = sd) {
  return (mean(x) - alpha*ineq_measure(x))
}

####............................................................................


#' Negative distance to a reference point
#'
#' @description
#' \code{negative_distance_to_reference_point} computes negative of distance in \eqn{L^p} norm between given vector and a reference point whose all coordinates are equal to a chosen aspiration level.
#'
#' @param x A numeric vector
#' @param A A real number
#' @param p A real number
#'
#' @details
#' Parameter \code{p} is the exponent defining the \eqn{L^p}-norm distance. The distance is properly defined if \eqn{p \geq 1}. If \eqn{p < 1} a warning message is displayed, but the value of the function is returned.
#'
#' Formula for negative distance to reference point:
#' \deqn{
#' g_{p, A}(x) = - \left( \sum_{i=1}^N |x_i - A|^p \right)^{\frac{1}{p}}
#' }
#' where \eqn{N} is the length of vector \eqn{x}.
#'
#' Parameter \eqn{A < \infty} is an aspiration level of utility. Parameter \eqn{p} defines the type of \eqn{L^p} distance and also determines the trade-off between equality (egaliatarianism for \eqn{p = \infty}) and efficiency (utilitarianism for \eqn{p = 1}).
#'
#' @return
#' Value of \eqn{g_{p, A}(x)}.
#' @export
#' negative_distance_to_reference_point
#'
#' @examples
#' negative_distance_to_reference_point(c(1,2,3,4))
#' negative_distance_to_reference_point(c(1,2,3,4), A = 100, p = 1.5)
#'
#'
negative_distance_to_reference_point <- function(x, A = 0, p  = 2) {
  if (p<1) {
    warning("L_p distance is not properly defined for p<1!")
  }

  return(-(sum(abs(x - A)^p))^(1/p))
}

####............................................................................


#' Log-sum-exp function
#'
#' @description
#' \code{log_sum_exp} function is a "soft min" of a vector (i.e., a differentiable function of a vector whose value is predominantly determined by its smallest element).
#'
#' @param x A numeric vector
#'
#' @details
#' Formula of log-sum-exp function:
#' \deqn{
#' g(x) = -\ln\left(\sum_{i=1}^{N} e^{-x_i} \right)
#' }
#'
#' @return
#' Value of \eqn{g(x)}.
#' @export
#' log_sum_exp
#'
#' @examples
#' log_sum_exp(c(1, 2, 3, 4))
#'
#'
log_sum_exp <- function(x) {
  return(-log(sum(exp(-x))))
}

####............................................................................


#' Minimum
#'
#' @description
#' \code{minimum} is the \eqn{\min(\cdot)} function.
#'
#' @param x A numeric vector
#'
#' @details
#' Formula:
#' \deqn{
#' g(x) = \min_{i = 1,\ldots,N} x_i
#' }
#' where \eqn{N} is the length of vector \code{x}.
#'
#' @return
#' Smallest element of the vector \code{x}.
#' @export
#' minimum
#'
#' @examples
#' minimum(c(1, 2, 3, 4))
#'
#'
minimum <- function(x) {
  return(min(x))
}

####............................................................................


#' Negative maximium downside semideviation
#'
#' @description
#' \code{negative_max_downside_semideviation} computes negative of the difference between mean of a vector and its smallest element.
#'
#' @param x A numeric vector
#'
#' @details
#' Formula:
#' \deqn{
#' g(x) = -\max_{i=1,\ldots,N} (\mu(x) - x_i)
#' }
#' where \eqn{N} is the length of vector \eqn{x} and \eqn{\mu(x) = \frac{1}{N}\sum_{i=1}^N x_i}.
#'
#' @return
#' Value of \eqn{g(x)}.
#' @export
#' negative_max_downside_semideviation
#'
#' @examples
#' negative_max_downside_semideviation(c(1,2,3,4))
#'
#'
negative_max_downside_semideviation <- function(x) {
  return(-max(mean(x) - x))
}

####............................................................................


#' Negative mean downside semideviation
#'
#' @description
#' \code{negative_mean_downside_semideviation} computes negative of the mean of differences between the mean of a vector and each of its elements that is smaller than the mean.
#'
#' @param x A numeric vector
#'
#' @details
#' Formula:
#' \deqn{
#' g(x) = -\frac{1}{N} \sum_{i=1}^N (\mu(x) - x_i)_+
#' }
#' where \eqn{N} is the length of vector \eqn{x}, \eqn{\mu(x) = \frac{1}{N}\sum_{i=1}^N x_i} and \eqn{(\cdot)_+ = \max\{\cdot, 0\}}.
#'
#' @return
#' Value of \eqn{g(x)}.
#' @export
#' negative_mean_downside_semideviation
#'
#' @examples
#' negative_mean_downside_semideviation(c(1,2,3,4))
#'
#'
negative_mean_downside_semideviation <- function(x) {
  N <- length(x)
  m <- mean(x)

  return(-sum((m - x)*(m - x > 0))/N)
}

####............................................................................


#' Negative total downside difference
#'
#' @description
#' \code{negative_total_downside_difference} computes negative of the sum of differences between each element of a vector and all other elements that are greater than it.
#'
#' @param x A numeric vector
#'
#' @details
#' Formula:
#' \deqn{
#' g(x) = - \sum_{i=1}^N \sum_{j=1}^N (x_j - x_i)_+
#' }
#' where \eqn{N} is the length of vector \eqn{x} and \eqn{(\cdot)_+ = \max\{\cdot, 0\}}.
#'
#' @return
#' Value of \eqn{g(x)}.
#' @export
#' negative_total_downside_difference
#'
#' @examples
#' negative_total_downside_difference(c(1,2,3,4))
#'
#'
negative_total_downside_difference <- function(x) {
  res <- 0

  for (u in x) {
    res <- res + sum((x-u)*(x-u > 0))
  }

  return(-res)
}

####............................................................................


#' Type-2 achievement function
#'
#' @description
#' \code{type_2_achievement_fun} computes a composition of minimum of elements of a vector (possibly in reference to a selected aspiration level) and the sum of elements of the vector.
#'
#' @param x A numeric vector
#' @param A A real number
#' @param alpha A real number
#'
#' @details
#' Parameter \code{alpha} should be non-negative, otherwise the function loses its interpretation. For negative values of \code{alpha} a warning message is displayed, but the value of the function is returned.
#'
#' Formula of Type-2 achievement function:
#' \deqn{
#' g_{\alpha, A}(x) = \min_{i=1,\ldots,N}(x_i - A) + \alpha \sum_{i=1}^N (x_i - A)
#' }
#' where \eqn{N} is the length of vector \eqn{x}.
#'
#' Type-2 achievement function is equivalent to a weighted average of vector \eqn{x} where its smallest element has a preferential weight of \eqn{(1 + 1/\alpha)} and all remaining elements have weight 1. Parameter \eqn{\alpha} determines the preference for equality (with \eqn{\alpha = 0} corresponding to strict egalitarianism and \eqn{\alpha = \infty} to utilitarianinsm).
#'
#' Parameter \eqn{A} represents an aspiration level of utility.
#'
#' @return
#' Value of \eqn{g_{\alpha, A}(x)}
#' @export
#' type_2_achievement_fun
#'
#' @examples
#' type_2_achievement_fun(c(1, 2, 3, 4))
#' type_2_achievement_fun(c(1, 2, 3, 4), A = 2, alpha = 0.2)
#'
#'
type_2_achievement_fun <- function(x, A = 0, alpha = 1) {
  if (alpha < 0) {
    warning("Type-2 achievement function has no interpretation when alpha < 0.")
  }
  return(min(x-A) + alpha*sum(x-A))
}

####............................................................................



#' Ordered weighted average (OWA)
#'
#' @description
#' OWA computes a weighted average of a vector that was sorted in an increasing order.
#'
#' @param x A numeric vector
#' @param w A numeric vector
#' @param normalize_weights Bool (TRUE by default)
#'
#' @details
#' All weights in the vector \code{w} should be non-negative for OWA to have a proper interpretation. In case of some weights being negative a warning message is displayed, but the value of the function is returned.
#'
#' Formula of OWA function:
#' \deqn{
#' g_w(x) = \sum_{i=1}^N w_i \ x_{(i)}
#' }
#' where \eqn{N} is the length of vector \eqn{x} and \eqn{x_{(i)}} denotes the \eqn{i}-th smallest element of \eqn{x}.
#'
#' Weights \eqn{w} forming a decreasing sequence of positive numbers, i.e., \eqn{w_1 > \ldots > w_N > 0} represent a decreasing priority to higher levels of utility in a utility profile \eqn{x}
#'
#' @return
#' Value of \eqn{g_w(x)}
#' @export
#' OWA
#'
#' @examples
#' OWA(c(1, 2, 3, 4), w = c(5, 4, 3, 2), normalize_weights = FALSE)
#'
#'
OWA <- function(x, w, normalize_weights = TRUE) {
  if (length(x) != length(w)) {
    stop("Vectors x and w must have equal lengths.")
  }

  if (any(w < 0)) {
    warning("OWA has no proper interpretation for some weithts being negative.")
  }
  u <- sort(x)
  z <- w

  if (normalize_weights == TRUE) {
    z <- z/sum(z)
  }

  return(u%*%z)
}

####............................................................................


#' Atkinson function
#'
#' @description
#' \code{Atkinson_fun} computes Atkinson function. It is also known as constant marginal utility of consumption function and has interpretation as prioritarian aggregating function.
#'
#' @param x A numeric vector
#' @param gamma A real number
#'
#' @details
#' Atkinson function is defined for any value of parameter \code{gamma}, however for negative values the function is convex and cannot be used as a prioritarian aggregating function.
#'
#' \code{x} is expected to be a utility profile, i.e., none of its elements is negative. If \code{x} contains negative values, a warning message is displayed, but the value of the function is returned.
#'
#' Formula of Atkinson function:
#' \deqn{
#' g_\gamma(x) =
#'  \left\{
#'    \begin{array}{ll}
#'      \sum_{i=1}^N \ln x_i &  \qquad \mathrm{for} \quad \gamma = 1\\
#'      (1-\gamma)^{-1} \sum_{i=1}^N x_i^{1-\gamma} &  \qquad \mathrm{for} \quad \gamma \neq 1
#'    \end{array}
#'  \right.
#' }
#' where \eqn{N} is the length of \eqn{x}.
#'
#' The parameter \eqn{\gamma} determines the social aversion to inequality: \eqn{\gamma = 0} corresponds to no aversion (utilitarianism), \eqn{\gamma = \infty} corresponds to full aversion to inequality (egalitarianism). The parameter \eqn{\gamma} can also be interpreted as socially tolerable degree of losses involved in transfers of utility from higher levels to lower levels.
#'
#' @return
#' Value of Atkinson function \eqn{g_\gamma(x)}.
#' @export
#' Atkinson_fun
#'
#' @examples
#' Atkinson_fun(c(1, 2, 3, 4))
#' Atkinson_fun(c(1, 2, 3, 4), gamma = 0.5)
#'
#'
Atkinson_fun <- function(x, gamma = 1) {
  if (any(x <= 0)) {
    warning("Atkinson function may be undefined for negative elements of x.")
  }

  res <- 0

  if (gamma == 1) {
    res <- sum(log(x))
  }
  else {
    res <- sum(x^(1-gamma))/(1-gamma)
  }

  return(res)
}

####............................................................................


#' Negative exponential function
#'
#' @description
#' \code{negative_exp_fun} computes negative of sum of negative exponents of elements of a vector.
#'
#' @param x A numeric vector
#'
#' @details
#' Formula of negative exponential function:
#' \deqn{
#' g(x) = -\sum_{i=1}^N e^{-x_i}
#' }
#' where \eqn{N} is the length of vector \eqn{x}.
#'
#' Negative exponential function is used as a prioritarian aggregating function that gives a very high priority to the smallest element of a utility vector \eqn{x}.
#'
#' @return
#' Value of \eqn{g(x)}.
#' @export
#' negative_exp_sum
#'
#' @examples
#' negative_exp_sum(c(1, 2, 3, 4))
#'
#'
negative_exp_sum <- function(x) {
  return(-sum(exp(-x)))
}


