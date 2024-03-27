#include "s21_math.h"

long double fact(double x);

long double s21_atan(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF)
    res = s21_PI / 2;
  else if (x == s21_NEG_INF)
    res = -s21_PI / 2;
  else {
    if (s21_fabs(x) < 1) {
      double tmp = x / s21_sqrt(1 + x * x);
      res = s21_asin(tmp);
    } else if (x > 1) {
      res = s21_PI / 2 - s21_atan(1 / x);
    } else {
      res = -s21_PI / 2 - s21_atan(1 / x);
    }
  }
  return res;
}

long double s21_fmod(double x, double y) {
  long double ans = 0, min_flag = 0;
  if ((x == S21_INF) || (x == S21_MIN_INF))
    ans = S21_NAN;
  else if (s21_fabs(y) < S21_EPS)
    ans = S21_NAN;
  else if ((y == S21_INF) || (y == S21_MIN_INF))
    ans = x;
  else {
    if (y < 0) {
      y *= -1;
    }
    if (x < 0) {
      x *= -1;
      min_flag = 1;
    }
    if (min_flag)
      ans = -1 * (x - y * (long long int)(x / y));
    else
      ans = x - y * (long long int)(x / y);
  }
  return ans;
}

long double s21_pow(double base, double exp) {
  long double ans = 1;
  if ((exp != exp) || (base != base)) {
    ans = S21_NAN;
  } else if (exp == -S21_INF) {
    if (base != 1)
      ans = 0;
    else
      ans = 1;
  } else if (exp == S21_INF) {
    if (base != 1)
      ans = S21_INF;
    else
      ans = 1;
  } else {
    if (exp > 0) {
      for (long long int i = 0; i < (long long int)exp; ++i) {
        ans *= base;
      }
    } else if (exp < 0) {
      for (long long int i = 0; i > (long long int)exp; --i) {
        ans /= base;
      }
    }
    long double expl = s21_fabs(exp) - s21_fabs((long long int)exp);
    if (expl > 0) {
      if (exp > 0) {
        ans *= s21_exp(expl * s21_log(base));
      } else if (exp < 0) {
        ans *= s21_exp(expl * s21_log(1 / base));
      }
    }
  }
  return ans;
}

long double s21_sqrt(double x) {
  long double S = x;
  long double xi = S / 2;
  long double xi_1 = S / 2;
  long double eps = 1;
  if (s21_fabs(S) < S21_EPS)
    xi = 0;
  else if (S < 0)
    xi = S21_NAN;
  else if (S != S)
    xi = S21_NAN;
  else if (S == S21_INF)
    xi = S21_INF;
  else if (S == S21_MIN_INF)
    xi = S21_NAN;
  else if (S > 0) {
    xi = S / 2;
    while ((S > 0) && (eps > S21_EPS_TAYLOR)) {
      xi = (xi_1 + S / xi_1) / 2;
      eps = -(S - xi * xi) / (2 * xi);
      xi_1 = xi;
    }
  }
  return xi;
}

long double s21_tan(double x) {
  long double res = 0;
  if (x == s21_POS_INF || x == s21_NEG_INF)
    res = -s21_nan;
  else if (x != x)
    res = s21_nan;
  else
    res = s21_sin(x) / s21_cos(x);
  return res;
}

long double s21_acos(double x) {
  long double res = 0;
  if (x != x || x > 1 || x < -1)
    res = s21_nan;
  else
    res = s21_PI / 2 - s21_asin(x);
  return res;
}

long double s21_asin(double x) {
  long double res = 0;
  if (x != x || x > 1 || x < -1)
    res = s21_nan;
  else {
    if (x == 1)
      res = s21_PI / 2;
    else if (x == -1)
      res = -s21_PI / 2;
    else {
      long double numerator = 1;
      long double denumerator = 2;
      long double elem = 1;
      long double poly = x * x * x;
      int current_degree = 3;
      while (s21_fabs(elem) >= epsilon) {
        elem = numerator * poly / denumerator / current_degree;
        res += elem;
        numerator *= current_degree;
        denumerator *= current_degree + 1;
        poly *= x * x;
        current_degree += 2;
      }
      res += x;
    }
  }
  return res;
}

long double s21_sin(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF || x == s21_NEG_INF)
    res = -s21_nan;
  else {
    x = s21_fmod(x, 100000 * s21_PI);
    long double numerator = x;
    long double denumerator = 1;
    long double elem = 1;
    int n = 1;
    while (s21_fabs(elem) >= epsilon) {
      elem = numerator / denumerator;
      res += elem;
      numerator *= -x * x;
      denumerator *= (n + 1) * (n + 2);
      n += 2;
    }
  }
  return res;
}

long double s21_cos(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF || x == s21_NEG_INF)
    res = -s21_nan;
  else {
    x = s21_fmod(x, 100000 * s21_PI);
    long double element = 0;
    long double stepen = x * x;
    long double factorial = 2;
    int count = 1;
    if (s21_fabs(x) <= 2 * s21_PI) {
      do {
        element = stepen / fact(factorial);
        if (count % 2 != 0) {
          res -= element;
        } else {
          res += element;
        }
        stepen *= x * x;
        factorial += 2;
        count++;
      } while (element >= epsilon);
      res += 1;
    } else if (x > 2 * s21_PI)
      res = s21_cos(x - 2 * s21_PI);
    else
      res = s21_cos(-x - 2 * s21_PI);
  }
  return res;
}

int is_int(double val) { return (val == ((double)((int)val))); }

long double s21_ceil(double x) {
  long double res;
  if (x != x) {
    res = S21_NAN;
  } else if (x == -S21_INF) {
    res = S21_INF;
  } else if (x == S21_INF) {
    res = S21_MIN_INF;
  } else if (is_int(x) || x == 0) {
    res = x;
  } else if (x > 0) {
    res = (int)x + 1;
  } else if (x < 0) {
    double dif = x - (int)x;
    res = x - dif;
  }
  return res;
}

long double s21_floor(double x) {
  long double ans = x >= 0. ? (int)x : (int)(x - 1.);
  if (x != x) {
    ans = S21_NAN;
  } else if (x == -S21_INF) {
    ans = S21_INF;
  } else if (x == S21_INF) {
    ans = S21_MIN_INF;
  } else if ((s21_fabs((int)x - x) < S21_EPS) && (x < 0))
    ++ans;
  return ans;
}

long double s21_fabs(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF || x == s21_NEG_INF)
    res = s21_POS_INF;
  else if (x < 0)
    res = -x;
  else
    res = x;
  return res;
}

int s21_abs(int x) {
  if (x < 0) x = -x;
  return x;
}

long double fact(double x) {
  long double result = 1;
  for (int i = x; i > 0; i--) {
    result *= i;
  }
  return result;
}

long double s21_exp(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x == s21_POS_INF)
    res = s21_POS_INF;
  else if (x == s21_NEG_INF)
    res = 0;
  else {
    res = 1;
    long double ans = 0;
    int n = 1;
    long double npow = x, nfac = 1;
    while (s21_fabs(res - ans) > 1e-15) {
      ans = res;
      res += npow / nfac;
      nfac *= ++n;
      npow *= x;
    }
  }
  return res;
}

long double s21_log(double x) {
  long double res = 0;
  if (x != x)
    res = s21_nan;
  else if (x < 0)
    res = -s21_nan;
  else if (x == s21_POS_INF)
    res = s21_POS_INF;
  else if (x == 0)
    res = s21_NEG_INF;
  else {
    long double temp = 1;
    long double xx = (x - 1) / (x + 1);
    long double stepen = xx;
    long double sum = 0;
    long double element = 0;
    do {
      element = stepen / temp;
      sum += element;
      temp += 2;
      stepen *= xx * xx;
    } while (s21_fabs(element) >= epsilon);
    res = 2 * sum;
  }
  return res;
}