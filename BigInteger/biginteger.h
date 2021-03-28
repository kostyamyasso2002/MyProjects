#include <iostream>
#include <vector>
#include <string>

class Rational;

class BigInteger {
 private:
  static const long long system_base = 100000;
  static const long long normal_system_base = 10;
  std::vector<long long> container_;
  bool is_positive_;

  void HelpNormalize(size_t border);
  void HardNormalize(size_t border);
  void Clear();
  void Increase(int x);
  void Decrease();
  void ChangeSign();

  friend class Rational;
  friend Rational operator-(const Rational& subtrahend);

 public:
  BigInteger();
  BigInteger(int x);
  BigInteger(const BigInteger& x);
  ~BigInteger() = default;


  BigInteger& operator++();
  BigInteger operator++(int);
  BigInteger& operator--();
  BigInteger operator--(int);
  BigInteger& operator+=(const BigInteger& added);
  BigInteger& operator-=(const BigInteger& added);
  friend BigInteger operator-(BigInteger& changed);

  BigInteger& operator*=(const BigInteger& multiplied);
  friend BigInteger operator*(const BigInteger& multiple_1, const BigInteger& multiple_2);

  BigInteger& operator/=(const BigInteger& divisor);
  BigInteger& operator%=(const BigInteger& divisor);

  friend bool operator<(const BigInteger& compared_1, const BigInteger& compared_2);

  BigInteger& operator=(BigInteger copy);
  void Swap(BigInteger& swapped);

  std::string toString() const;
  friend std::ostream& operator<<(std::ostream& out, const BigInteger& printed);
  friend std::istream& operator>>(std::istream& in, BigInteger& input);

  explicit operator bool() const;
  explicit operator double() const;
};


class Rational {
 private:
  BigInteger numerator_;
  BigInteger denominator_;

  BigInteger Gcd(BigInteger x, BigInteger y);
  void NormalizeSign();
  void Normalize();
  BigInteger FastPow(BigInteger x, long long st) const;
 public:
  Rational();
  Rational(int x);
  Rational(const BigInteger& x);
  Rational(const BigInteger& numerator, const BigInteger& denominator);
  ~Rational() = default;

  Rational& operator+=(const Rational& term);

  Rational& operator-=(const Rational& subtrahend);
  friend Rational operator-(const Rational& subtrahend);

  Rational& operator*=(const Rational& multiply);

  Rational& operator/=(const Rational& divisor);

  friend bool operator<(const Rational& lhs, const Rational& rhs);

  std::string toString() const;
  std::string asDecimal(size_t precision = 0) const;

  explicit operator double() const;
};

void BigInteger::HelpNormalize(size_t border) {
  for (size_t i = 0; i < container_.size(); ++i) {
    if (i == border && border != 0 && container_[i] >=0 && container_[i] < system_base) {
      return;
    }
    if (container_[i] < 0) {
      if (i == container_.size() - 1) {
        continue;
      }
      long long t = (-container_[i]) / system_base + 1;
      container_[i] += t * system_base;
      container_[i + 1] -= t;
    }
    if (container_[i] >= system_base) {
      if (i == container_.size() - 1) {
        container_.push_back(container_[i] / system_base);
      } else {
        container_[i + 1] += container_[i] / system_base;
      }
      container_[i] %= system_base;
    }
  }

  for (size_t i = container_.size() - 1; i >= 1; --i) {
    if (container_[i] == 0) {
      container_.pop_back();
    } else {
      break;
    }
  }
}

void BigInteger::HardNormalize(size_t border = 0) {
  HelpNormalize(border);

  if (container_.back() < 0) {
    is_positive_ = !is_positive_;
    for (size_t i = 0; i < container_.size(); ++i) {
      container_[i] *= -1;
    }
    HelpNormalize(border);
  }

  if (container_.size() == 1 && container_[0] == 0) {
    is_positive_ = true;
  }
}

void BigInteger::Clear() {
  is_positive_ = true;
  container_.clear();
}

void BigInteger::ChangeSign() {
  is_positive_ = !is_positive_;
  if (container_.size() == 1 && container_[0] == 0) {
    is_positive_ = true;
  }
}

void BigInteger::Increase(int x) {
  container_.resize(container_.size() + x);
  for (int i = container_.size() - 1; i >= 0; --i) {
    if (i >= x) {
      container_[i] = container_[i - x];
    } else {
      container_[i] = 0;
    }
  }
  HardNormalize();
}

void BigInteger::Decrease() {
  for (size_t i = 0; i < container_.size() - 1; ++i) {
    container_[i] = container_[i + 1];
  }
  container_[container_.size() - 1] = 0;
  HardNormalize();
}

BigInteger::BigInteger(): BigInteger(0) {}

BigInteger::BigInteger(int x) {
  if (x < 0) {
    x *= -1;
    is_positive_ = false;
  } else {
    is_positive_ = true;
  }
  while (x > 0) {
    container_.push_back(x % system_base);
    x /= system_base;
  }
  if (container_.empty()) {
    container_.push_back(0);
  }
}

BigInteger::BigInteger(const BigInteger& x) {
  container_ = x.container_;
  is_positive_ = x.is_positive_;
}

BigInteger& BigInteger::operator++() {
  if (is_positive_) {
    ++container_[0];
    HardNormalize(1);
    return *this;
  }
  --container_[0];
  HardNormalize(1);
  return *this;
}

BigInteger BigInteger::operator++(int) {
  BigInteger copy = *this;
  ++*this;
  return copy;
}

BigInteger& BigInteger::operator--() {
  if (is_positive_) {
    --container_[0];
    HardNormalize(1);
    return *this;
  }
  ++container_[0];
  HardNormalize(1);
  return *this;
}

BigInteger BigInteger::operator--(int) {
  BigInteger copy = *this;
  --*this;
  return copy;
}

BigInteger& BigInteger::operator+=(const BigInteger& added) {
  if (container_.size() < added.container_.size()) {
    container_.resize(added.container_.size());
  }
  for (size_t i = 0; i < added.container_.size(); ++i) {
    if (added.is_positive_ == is_positive_) {
      container_[i] += added.container_[i];
    } else {
      container_[i] -= added.container_[i];
    }
  }
  HardNormalize(added.container_.size() + 1);
  return *this;
}


BigInteger& BigInteger::operator-=(const BigInteger& added) {
  if (container_.size() < added.container_.size()) {
    container_.resize(added.container_.size());
  }
  for (size_t i = 0; i < added.container_.size(); ++i) {
    if (added.is_positive_ == is_positive_) {
      container_[i] -= added.container_[i];
    } else {
      container_[i] += added.container_[i];
    }
  }
  HardNormalize(added.container_.size() + 1);
  return *this;
}

BigInteger operator+(const BigInteger& term_1, const BigInteger& term_2) {
  BigInteger copy = term_1;
  copy += term_2;
  return copy;
}

BigInteger operator-(const BigInteger& term_1, const BigInteger& term_2) {
  BigInteger copy = term_1;
  copy -= term_2;
  return copy;
}

BigInteger operator-(BigInteger& changed) {
  BigInteger copy = changed;
  copy.ChangeSign();
  return copy;
}

BigInteger& BigInteger::operator*=(const BigInteger& multiplied) {
  is_positive_ = is_positive_ == multiplied.is_positive_;
  container_.resize(container_.size() + multiplied.container_.size());
  for (size_t i = container_.size() - 1; true; --i) {
    long long value = 0;
    for (size_t j = 0; j < multiplied.container_.size() && i >= j; ++j) {
      value += multiplied.container_[j] * container_[i - j];
    }
    container_[i] = value;
    if (i == 0) {
      break;
    }
  }
  HardNormalize();
  return *this;
}

BigInteger operator*(const BigInteger& multiple_1, const BigInteger& multiple_2) {
  BigInteger copy = multiple_1;
  copy *= multiple_2;
  return copy;
}

BigInteger& BigInteger::operator/=(const BigInteger& divisor) {
  BigInteger result;
  result.is_positive_ = (is_positive_ == divisor.is_positive_);
  BigInteger divisor_copy = divisor;
  divisor_copy.is_positive_ = true;
  is_positive_ = true;

  result.container_.resize(container_.size());
  if (container_.size() < divisor.container_.size()) {
    *this = BigInteger(0);
    return *this;
  }
  divisor_copy.Increase(container_.size() - divisor.container_.size());
  for (int i = container_.size() - divisor.container_.size(); i >= 0; --i) {
    long long left = 0;
    long long right = system_base;
    while (right - left > 1) {
      long long mid = (left + right) / 2;
      if (!(*this < divisor_copy * mid)) {
        left = mid;
      } else {
        right = mid;
      }
    }
    *this -= left * divisor_copy;
    result.container_[i] += left;
    divisor_copy.Decrease();
  }
  result.HardNormalize();
  *this = result;
  HardNormalize();
  return *this;
}

BigInteger& BigInteger::operator%=(const BigInteger& divisor) {
  bool result_sign = is_positive_;
  BigInteger divisor_copy = divisor;
  divisor_copy.is_positive_ = true;
  is_positive_ = true;

  if (container_.size() < divisor.container_.size()) {
    is_positive_ = result_sign;
    return *this;
  }
  divisor_copy.Increase(container_.size() - divisor.container_.size());
  for (int i = container_.size() - divisor.container_.size(); i >= 0; --i) {
    long long left = 0;
    long long right = system_base;
    while (right - left > 1) {
      long long mid = (left + right) / 2;
      if (!(*this < divisor_copy * mid)) {
        left = mid;
      } else {
        right = mid;
      }
    }
    *this -= left * divisor_copy;
    divisor_copy.Decrease();
  }

  is_positive_ = result_sign;
  HardNormalize();
  return *this;
}

BigInteger operator/(const BigInteger& divided, const BigInteger& divisor) {
  BigInteger copy = divided;
  copy /= divisor;
  return copy;
}


BigInteger operator%(const BigInteger& divided, const BigInteger& divisor) {
  BigInteger copy = divided;
  copy %= divisor;
  return copy;
}

bool operator<(const BigInteger& compared_1, const BigInteger& compared_2) {
  if (compared_1.is_positive_ != compared_2.is_positive_) {
    return !compared_1.is_positive_;
  }
  if (compared_1.container_.size() != compared_2.container_.size()) {
    return compared_1.is_positive_ == (compared_1.container_.size() < compared_2.container_.size());
  }
  for (int i = compared_2.container_.size() - 1; i >= 0; --i) {
    if (compared_1.container_[i] != compared_2.container_[i]) {
      return compared_1.is_positive_ == (compared_1.container_[i] < compared_2.container_[i]);
    }
  }
  return false;
}

bool operator>(const BigInteger& compared_1, const BigInteger& compared_2) {
  return compared_2 < compared_1;
}

bool operator>=(const BigInteger& compared_1, const BigInteger& compared_2) {
  return !(compared_1 < compared_2);
}

bool operator<=(const BigInteger& compared_1, const BigInteger& compared_2) {
  return !(compared_2 < compared_1);
}

bool operator==(const BigInteger& compared_1, const BigInteger& compared_2) {
  return !(compared_1 < compared_2 || compared_2 < compared_1);
}

bool operator!=(const BigInteger& compared_1, const BigInteger& compared_2) {
  return compared_1 < compared_2 || compared_2 < compared_1;
}

BigInteger& BigInteger::operator=(BigInteger copy) {
  Swap(copy);
  return *this;
}

void BigInteger::Swap(BigInteger& swapped) {
  std::swap(is_positive_, swapped.is_positive_);
  std::swap(container_, swapped.container_);
}

std::string BigInteger::toString() const {
  std::string returned;
  if (!is_positive_) {
    returned += "-";
  }
  for (size_t i = container_.size() - 1; true; --i) {
    if (i == container_.size() - 1) {
      returned += std::to_string(container_[i]);
      if (i == 0) {
        break;
      }
      continue;
    }
    long long k = 0;
    long long st = 1;
    while (std::max(container_[i], 1LL) * st < BigInteger::system_base) {
      ++k;
      st *= BigInteger::normal_system_base;
    }
    for (int j = 0; j < k - 1; ++j) {
      returned += "0";
    }
    returned += std::to_string(container_[i]);
    if (i == 0) {
      break;
    }
  }
  return returned;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& printed) {
  out << printed.toString();
  return out;
}

std::istream& operator>>(std::istream& in, BigInteger& input) {
  std::string input_string;
  in >> input_string;
  input.Clear();
  if (input_string[0] == '-') {
    input.is_positive_ = false;
  }
  long long multiple = 1;
  long long digit = 0;
  input.container_.resize(1);
  for (size_t i = input_string.size() - 1; true; --i) {
    if (input_string[i] == '-' || input_string[i] == '+') {
      break;
    }
    if (multiple == BigInteger::system_base) {
      multiple = 1;
      ++digit;
      input.container_.push_back(0);
    }
    input.container_[digit] += multiple * (input_string[i] - '0');
    multiple *= BigInteger::normal_system_base;
    if (i == 0) {
      break;
    }
  }
  input.HardNormalize();
  return in;
}

BigInteger::operator bool() const {
  return *this != BigInteger(0);
}

BigInteger::operator double() const {
  double help = 1;
  double result = 0;
  for (size_t i = 0; i < container_.size(); ++i) {
    result += help * container_[i];
    help *= system_base;
  }
  return result;
}


BigInteger operator "" _bi(unsigned long long x) {
  return BigInteger(x);
}

BigInteger Rational::Gcd(BigInteger x, BigInteger y) {
  x.is_positive_ = true;
  y.is_positive_ = true;
  if (x == 0) {
    return y;
  }
  if (y == 0) {
    return x;
  }
  if (x < y) {
    return Gcd(x, y % x);
  }
  return Gcd(y, x % y);
}

void Rational::NormalizeSign() {
  if (!denominator_.is_positive_) {
    denominator_.is_positive_ = true;
    numerator_.ChangeSign();
  }
}

void Rational::Normalize() {
  NormalizeSign();
  BigInteger t = Gcd(numerator_, denominator_);
  numerator_ /= t;
  denominator_ /= t;
}

BigInteger Rational::FastPow(BigInteger x, long long st) const {
  if (st == 0) {
    return BigInteger(1);
  }
  BigInteger res = 1;
  if (st % 2 == 1) {
    res = x;
  }
  BigInteger ret = FastPow(x, st / 2);
  res *= ret;
  res *= ret;
  return res;
}

Rational::Rational(): numerator_(BigInteger(0)), denominator_(BigInteger(1)) {}

Rational::Rational(int x): numerator_(BigInteger(x)), denominator_(BigInteger(1)) {}

Rational::Rational(const BigInteger& x): numerator_(x), denominator_(BigInteger(1)) {}

Rational::Rational(const BigInteger& numerator, const BigInteger& denominator):
        numerator_(numerator), denominator_(denominator) {}

Rational& Rational::operator+=(const Rational& term) {
  numerator_ *= term.denominator_;
  numerator_ += denominator_ * term.numerator_;
  denominator_ *= term.denominator_;
  Normalize();
  return *this;
}

Rational operator+(const Rational& term_1, const Rational& term_2) {
  Rational copy = term_1;
  copy += term_2;
  return copy;
}

Rational& Rational::operator-=(const Rational& subtrahend) {
  numerator_ *= subtrahend.denominator_;
  numerator_ -= denominator_ * subtrahend.numerator_;
  denominator_ *= subtrahend.denominator_;
  Normalize();
  return *this;
}

Rational operator-(const Rational& subtrahend_1, const Rational& subtrahend_2) {
  Rational copy = subtrahend_1;
  copy -= subtrahend_2;
  return copy;
}

Rational operator-(const Rational& subtrahend) {
  Rational copy = subtrahend;
  copy.numerator_.ChangeSign();
  return copy;
}

Rational& Rational::operator*=(const Rational& multiply) {
  numerator_ *= multiply.numerator_;
  denominator_ *= multiply.denominator_;
  Normalize();
  return *this;
}

Rational operator*(const Rational& multiply_1, const Rational& multiply_2) {
  Rational copy = multiply_1;
  copy *= multiply_2;
  return copy;
}

Rational& Rational::operator/=(const Rational& divisor) {
  numerator_ *= divisor.denominator_;
  denominator_ *= divisor.numerator_;
  Normalize();
  return *this;
}

Rational operator/(const Rational& divided, const Rational& divisor) {
  Rational copy = divided;
  copy /= divisor;
  return copy;
}

bool operator<(const Rational& lhs, const Rational& rhs) {
  return lhs.numerator_ * rhs.denominator_ < lhs.denominator_ * rhs.numerator_;
}

bool operator>(const Rational& lhs, const Rational& rhs) {
  return rhs < lhs;
}

bool operator>=(const Rational& lhs, const Rational& rhs) {
  return !(lhs < rhs);
}

bool operator<=(const Rational& lhs, const Rational& rhs) {
  return !(rhs < lhs);
}

bool operator==(const Rational& lhs, const Rational& rhs) {
  return !(lhs < rhs || rhs < lhs);
}

bool operator!=(const Rational& lhs, const Rational& rhs) {
  return (lhs < rhs || rhs < lhs);
}

std::string Rational::toString() const {
  std::string returned;
  returned += numerator_.toString();
  if (denominator_ != 1) {
    returned += '/';
    returned += denominator_.toString();
  }
  return returned;
}

std::string Rational::asDecimal(size_t precision) const {
  BigInteger st = 1;
  st = FastPow(BigInteger::normal_system_base, precision);
  std::string returned;
  BigInteger result = numerator_ * st / denominator_;
  if (numerator_ >= 0) {
    returned = result.toString();
  } else {
    returned = (-result).toString();
  }
  if (precision == 0) {
    if (numerator_ < 0) {
      returned = "-" + returned;
    }
    return returned;
  }
  while (returned.size() <= precision) {
    returned = "0" + returned;
  }
  returned = "0" + returned;
  for (size_t i = 0; i < returned.size() - precision - 1; ++i) {
    returned[i] = returned[i + 1];
  }
  returned[returned.size() - precision - 1] = '.';

  if (numerator_ < 0) {
    returned = "-" + returned;
  }
  return returned;
}

Rational::operator double() const {
  double num = static_cast<double>(numerator_);
  double den = static_cast<double>(denominator_);

  return num / den;
}
