#include <vector>
#include <cmath>

bool Equal(double a, double b) {
  return (std::abs(a - b) < 0.0000001);
}

class Line;

struct Point {
 public:
  double x;
  double y;

  Point(double x, double y);
  Point();
  ~Point() = default;

  friend bool operator==(const Point& first, const Point& second);
  friend bool operator!=(const Point& first, const Point& second);

  Point rotate(const Point& center, double angle) const;
  Point reflex(const Point& center) const;
  Point reflex(const Line& axis) const;
  Point scale(const Point& center, double coefficient) const;
};

class Line {
 private:
  //прямая Ax + By + C = 0
  double A;
  double B;
  double C;

 public:
  Line(const Point& first, const Point& second);
  Line(double k, double b);
  Line(const Point& point, double k);
  ~Line() = default;

  friend bool operator==(const Line& lhs, const Line& rhs);
  friend bool operator!=(const Line& lhs, const Line& rhs);
  friend Point Point::reflex(const Line& axis) const;
  friend Point Intersection(const Line& line_1, const Line& line_2);
};

class FlatVector {
 private:
  double x;
  double y;

 public:
  FlatVector(const Point& first, const Point& second);
  ~FlatVector() = default;

  friend double Scalar(const FlatVector& first, const FlatVector& second);
  friend double Cross(const FlatVector& first, const FlatVector& second);

  double Length() const;
};

class Shape {
 public:
  virtual ~Shape() = default;

  virtual double perimeter() const = 0;
  virtual double area() const = 0;

  virtual bool operator==(const Shape& another) const = 0;
  virtual bool operator!=(const Shape& another) const = 0;

  virtual bool isCongruentTo(const Shape& another) const = 0;
  virtual bool isSimilarTo(const Shape& another) const = 0;
  virtual bool containsPoint(const Point& point) const = 0;

  virtual void rotate(const Point& center, double angle) = 0;
  virtual void reflex(const Point& center) = 0;
  virtual void reflex(const Line& axis) = 0;
  virtual void scale(const Point& center, double coefficient) = 0;
};

class Polygon : public Shape {
 protected:
  std::vector<Point> points_;

  static bool AreEqualStraight(const std::vector<std::pair<double, double>>& first,
                               const std::vector<std::pair<double, double>>& second);
  static bool AreEqualBack(const std::vector<std::pair<double, double>>& first,
                           const std::vector<std::pair<double, double>>& second);
  Point Center() const;
  std::vector<std::pair<double, double>> LengthsAndAnglesRight() const;
  std::vector<std::pair<double, double>> LengthsAndAnglesLeft() const;

 public:
  Polygon(const std::vector<Point>& points);
  Polygon(std::initializer_list<Point> points_list);
  virtual ~Polygon() override = default;

  int verticesCount();
  const std::vector<Point>& getVertices();
  bool isConvex();

  virtual double perimeter() const override;
  virtual double area() const override;

  virtual bool operator==(const Shape& another) const override;
  virtual bool operator!=(const Shape& another) const override;

  virtual bool isCongruentTo(const Shape& another) const override;
  virtual bool isSimilarTo(const Shape& another) const override;
  virtual bool containsPoint(const Point& point) const override;

  virtual void rotate(const Point& center, double angle) override;
  virtual void reflex(const Point& center) override;
  virtual void reflex(const Line& axis) override;
  virtual void scale(const Point& center, double coefficient) override;
};

class Ellipse : public Shape {
 private:
  Point focus_1_;
  Point focus_2_;

 protected:
  double distance_;
  double small_semi_axis_;
  double big_semi_axis_;

 public:
  Ellipse(const Point& focus_1, const Point& focus_2, double double_distance);
  virtual ~Ellipse() override = default;

  std::pair<Point, Point> focuses() const;
  std::pair<Line, Line> directrices() const;
  double eccentricity() const;
  Point center() const;

  virtual double perimeter() const override;
  virtual double area() const override;

  virtual bool operator==(const Shape& another) const override;
  virtual bool operator!=(const Shape& another) const override;

  virtual bool isCongruentTo(const Shape& another) const override;
  virtual bool isSimilarTo(const Shape& another) const override;
  virtual bool containsPoint(const Point& point) const override;

  virtual void rotate(const Point& center, double angle) override;
  virtual void reflex(const Point& center) override;
  virtual void reflex(const Line& axis) override;
  virtual void scale(const Point& center, double coefficient) override;
};

class Circle : public Ellipse {
 public:
  Circle(const Point& center, double radius);
  virtual ~Circle() override = default;

  double radius() const;
};

class Rectangle : public Polygon {
 public:
  Rectangle(const Point& first, const Point& second, double ratio);
  virtual ~Rectangle() override = default;

  Point center();
  std::pair<Line, Line> diagonals();
};

class Square : public Rectangle {
 public:
  Square(const Point& first, const Point& second);
  virtual ~Square() override = default;

  Circle circumscribedCircle();
  Circle inscribedCircle();
};

class Triangle : public Polygon {
 public:
  Triangle(const Point& first, const Point& second, const Point& third);
  virtual ~Triangle() override = default;

  Circle circumscribedCircle();
  Circle inscribedCircle();

  Point centroid();
  Point orthocenter();
  Line EulerLine();
  Circle ninePointsCircle();
};

Point::Point(double x, double y): x(x), y(y) {}

Point::Point(): Point(0, 0) {}

bool operator==(const Point& first, const Point& second) {
  return (Equal(first.x, second.x)) && (Equal(first.y, second.y));
}

bool operator!=(const Point& first, const Point& second) {
  return !(first == second);
}

Point Point::rotate(const Point& center, double angle) const {
  angle *= M_PI / 180;
  return Point((x - center.x) * cos(angle) - (y - center.y) * sin(angle) + center.x,
               (x - center.x) * sin(angle) + (y - center.y) * cos(angle) + center.y);
}

Point Point::reflex(const Point& center) const {
  return Point(2 * center.x - x, 2 * center.y - y);
}

Point Point::reflex(const Line& axis) const {
  return Point((x * (axis.B * axis.B - axis.A * axis.A) - 2 * y * axis.A * axis.B -
                2 * axis.A * axis.C) / (axis.A * axis.A + axis.B * axis.B),
               (y * (axis.A * axis.A - axis.B * axis.B) - 2 * x * axis.A * axis.B -
                2 * axis.B * axis.C) / (axis.A * axis.A + axis.B * axis.B));
}

Point Intersection(const Line& line_1, const Line& line_2) {
  return Point(
          (line_2.C * line_1.B - line_2.B * line_1.C) / (line_1.A * line_2.B - line_2.A * line_1.B),
          (line_1.A * line_2.C - line_1.C * line_2.A) /
          (line_1.B * line_2.A - line_1.A * line_2.B));
}


Point Point::scale(const Point& center, double coefficient) const {
  return Point((x - center.x) * coefficient + center.x, (y - center.y) * coefficient + center.y);
}

FlatVector::FlatVector(const Point& first, const Point& second): x(second.x - first.x),
                                                                 y(second.y - first.y) {}

double Scalar(const FlatVector& first, const FlatVector& second) {
  return first.x * second.x + first.y * second.y;
}

double Cross(const FlatVector& first, const FlatVector& second) {
  return first.x * second.y - first.y * second.x;
}

double FlatVector::Length() const {
  return sqrt(x * x + y * y);
}


Line::Line(const Point& first, const Point& second): A(first.y - second.y), B(second.x - first.x),
                                                     C(first.x * second.y - second.x * first.y) {}


Line::Line(double k, double b): A(k), B(-1), C(b) {}


Line::Line(const Point& point, double k): Line(k, point.y - k * point.x) {}

bool operator==(const Line& lhs, const Line& rhs) {
  if (!Equal(lhs.A, 0) || !Equal(rhs.A, 0)) {
    if (Equal(lhs.A, 0) || Equal(rhs.A, 0)) {
      return false;
    }
    return (Equal(lhs.B * rhs.A, rhs.B * lhs.A)) && (Equal(lhs.C * rhs.A, rhs.C * lhs.A));
  }
  if (Equal(lhs.B, 0) || Equal(rhs.B, 0)) {
    return false;
  }
  return Equal(lhs.C * rhs.B, rhs.C * lhs.B);
}

bool operator!=(const Line& lhs, const Line& rhs) {
  return !(lhs == rhs);
}


bool Polygon::AreEqualStraight(const std::vector<std::pair<double, double>>& first,
                               const std::vector<std::pair<double, double>>& second) {
  if (first.size() != second.size()) {
    return false;
  }

  for (size_t i = 0; i < first.size(); ++i) {
    bool is_equal = true;
    for (size_t j = 0; j < first.size(); ++j) {
      if (!Equal(first[j].first, second[(i + j) % first.size()].first) ||
          !Equal(first[j].second, second[(i + j) % first.size()].second)) {
        is_equal = false;
        break;
      }
    }
    if (is_equal) {
      return true;
    }
  }

  return false;
}

bool Polygon::AreEqualBack(const std::vector<std::pair<double, double>>& first,
                           const std::vector<std::pair<double, double>>& second) {
  if (first.size() != second.size()) {
    return false;
  }

  for (size_t i = 0; i < first.size(); ++i) {
    bool is_equal = true;
    for (size_t j = 0; j < first.size(); ++j) {
      if (!Equal(first[j].first, second[(i + first.size() - j) % first.size()].first) ||
          !Equal(first[j].second, second[(i + first.size() - j) % first.size()].second)) {
        is_equal = false;
        break;
      }
    }
    if (is_equal) {
      return true;
    }
  }

  return false;
}

Point Polygon::Center() const {
  double x = 0;
  double y = 0;
  for (size_t i = 0; i < points_.size(); ++i) {
    x += points_[i].x;
    y += points_[i].y;
  }
  return Point(x / points_.size(), y / points_.size());
}

std::vector<std::pair<double, double>> Polygon::LengthsAndAnglesLeft() const {
  std::vector<std::pair<double, double>> answer;
  for (size_t i = 0; i < points_.size(); ++i) {
    FlatVector side_1(points_[i], points_[(i + points_.size() - 1) % points_.size()]);
    FlatVector side_2(points_[i], points_[(i + 1) % points_.size()]);
    answer.push_back({side_1.Length(), Scalar(side_1, side_2)});
  }
  return answer;
}

std::vector<std::pair<double, double>> Polygon::LengthsAndAnglesRight() const {
  std::vector<std::pair<double, double>> answer;
  for (size_t i = 0; i < points_.size(); ++i) {
    FlatVector side_1(points_[i], points_[(i + points_.size() - 1) % points_.size()]);
    FlatVector side_2(points_[i], points_[(i + 1) % points_.size()]);
    answer.push_back({side_2.Length(), Scalar(side_1, side_2)});
  }
  return answer;
}

Polygon::Polygon(const std::vector<Point>& points): points_(points) {}

Polygon::Polygon(std::initializer_list<Point> points_list): points_(points_list) {}

int Polygon::verticesCount() {
  return points_.size();
}

const std::vector<Point>& Polygon::getVertices() {
  return points_;
}

bool Polygon::isConvex() {
  int small_angle = 0;
  int big_angle = 0;
  for (size_t i = 0; i < points_.size(); ++i) {
    Point point_1 = points_[(i + points_.size() - 1) % points_.size()];
    Point point_2 = points_[i];
    Point point_3 = points_[(i + 1) % points_.size()];
    FlatVector side_1(point_2, point_1);
    FlatVector side_2(point_2, point_3);
    if (Cross(side_1, side_2) > 0) {
      ++small_angle;
    } else {
      ++big_angle;
    }
  }
  return small_angle == 0 || big_angle == 0;
}

double Polygon::perimeter() const {
  double answer = 0;
  for (size_t i = 0; i < points_.size(); ++i) {
    Point point_1 = points_[i];
    Point point_2 = points_[(i + 1) % points_.size()];
    FlatVector side(point_1, point_2);
    answer += side.Length();
  }
  return answer;
}

double Polygon::area() const {
  double answer = 0;
  for (size_t i = 0; i < points_.size(); ++i) {
    Point point_1 = points_[i];
    Point point_2 = points_[(i + 1) % points_.size()];
    answer += point_1.x * point_2.y - point_2.x * point_1.y;
  }
  return std::abs(answer) / 2;
}

bool Polygon::operator==(const Shape& another) const {
  if (dynamic_cast<const Polygon*>(&another) == nullptr) {
    return false;
  }
  const Polygon& other = dynamic_cast<const Polygon&>(another);

  std::vector<std::pair<double, double>> points_1;
  std::vector<std::pair<double, double>> points_2;
  for (size_t i = 0; i < points_.size(); ++i) {
    points_1.push_back({points_[i].x, points_[i].y});
  }
  for (size_t i = 0; i < other.points_.size(); ++i) {
    points_2.push_back({other.points_[i].x, other.points_[i].y});
  }

  return AreEqualStraight(points_1, points_2) || AreEqualBack(points_1, points_2);
}

bool Polygon::operator!=(const Shape& another) const {
  return !(*this == another);
}

bool Polygon::isCongruentTo(const Shape& another) const {
  if (dynamic_cast<const Polygon*>(&another) == nullptr) {
    return false;
  }
  const Polygon& other = dynamic_cast<const Polygon&>(another);
  return (AreEqualStraight(LengthsAndAnglesLeft(), other.LengthsAndAnglesLeft()) ||
          AreEqualBack(LengthsAndAnglesLeft(), other.LengthsAndAnglesRight()));
}

bool Polygon::isSimilarTo(const Shape& another) const {
  if (dynamic_cast<const Polygon*>(&another) == nullptr) {
    return false;
  }
  const Polygon& other = dynamic_cast<const Polygon&>(another);

  std::vector<std::pair<double, double>> lengths_1 = LengthsAndAnglesRight();
  std::vector<std::pair<double, double>> lengths_2 = other.LengthsAndAnglesRight();
  std::vector<std::pair<double, double>> lengths_3 = LengthsAndAnglesLeft();

  double perimeter_1 = perimeter();
  double perimeter_2 = other.perimeter();
  for (size_t i = 0; i < points_.size(); ++i) {
    lengths_1[i].first /= perimeter_1;
    lengths_1[i].second /= perimeter_1 * perimeter_1;
    lengths_3[i].first /= perimeter_1;
    lengths_3[i].second /= perimeter_1 * perimeter_1;
  }
  for (size_t i = 0; i < other.points_.size(); ++i) {
    lengths_2[i].first /= perimeter_2;
    lengths_2[i].second /= perimeter_2 * perimeter_2;
  }
  return (AreEqualStraight(lengths_1, lengths_2) || AreEqualBack(lengths_2, lengths_3));
}

bool Polygon::containsPoint(const Point& point) const {
  bool answer = false;
  for (size_t i = 0; i < points_.size(); ++i) {
    Point point_1 = points_[i].rotate(point, 1);
    Point point_2 = points_[(i + 1) % points_.size()].rotate(point, 1);
    if (((point_1.y < point.y && point_2.y >= point.y) ||
         (point_2.y < point.y && point_1.y >= point.y)) &&
        (point_1.x + (point.y - point_1.y) / (point_2.y - point_1.y) * (point_2.x - point_1.x) <
         point.x)) {
      answer = !answer;
    }
  }
  return answer;
}

void Polygon::rotate(const Point& center, double angle) {
  for (size_t i = 0; i < points_.size(); ++i) {
    points_[i] = points_[i].rotate(center, angle);
  }
}

void Polygon::reflex(const Point& center) {
  for (size_t i = 0; i < points_.size(); ++i) {
    points_[i] = points_[i].reflex(center);
  }
}

void Polygon::reflex(const Line& axis) {
  for (size_t i = 0; i < points_.size(); ++i) {
    points_[i] = points_[i].reflex(axis);
  }
}

void Polygon::scale(const Point& center, double coefficient) {
  for (size_t i = 0; i < points_.size(); ++i) {
    points_[i] = points_[i].scale(center, coefficient);
  }
}

Ellipse::Ellipse(const Point& focus_1, const Point& focus_2, double double_distance):
        focus_1_(focus_1), focus_2_(focus_2), distance_(double_distance / 2) {
  small_semi_axis_ = sqrt(
          distance_ * distance_ - pow(FlatVector(focus_1_, focus_2_).Length() / 2, 2));
  big_semi_axis_ = distance_;
}

std::pair<Point, Point> Ellipse::focuses() const {
  return {focus_1_, focus_2_};
}

std::pair<Line, Line> Ellipse::directrices() const {
  double coefficient = eccentricity() * eccentricity();
  Point point_1 = focus_1_.scale(center(), 1 / coefficient);
  Point point_2(point_1.x + (focus_1_.y - focus_2_.y), point_1.y - (focus_1_.x - focus_2_.x));
  Line dir_1(point_1, point_2);
  Line dir_2(point_1.reflex(center()), point_2.reflex(center()));
  return {dir_1, dir_2};
}

double Ellipse::eccentricity() const {
  return FlatVector(focus_1_, focus_2_).Length() / (2 * distance_);
}

Point Ellipse::center() const {
  return Point((focus_1_.x + focus_2_.x) / 2, (focus_1_.y + focus_2_.y) / 2);
}

double Ellipse::perimeter() const {
  return 4 * big_semi_axis_ * std::comp_ellint_2(eccentricity());
}

double Ellipse::area() const {
  return M_PI * small_semi_axis_ * big_semi_axis_;
}


bool Ellipse::operator==(const Shape& another) const {
  if (dynamic_cast<const Ellipse*>(&another) == nullptr) {
    return false;
  }
  const Ellipse& other = dynamic_cast<const Ellipse&>(another);
  if (!Equal(distance_, other.distance_)) {
    return false;
  }
  return ((focus_1_ == other.focus_2_ && focus_2_ == other.focus_1_) ||
          (focus_1_ == other.focus_1_ && focus_2_ == other.focus_2_));
}

bool Ellipse::operator!=(const Shape& another) const {
  return !(*this == another);
}

bool Ellipse::isCongruentTo(const Shape& another) const {
  if (dynamic_cast<const Ellipse*>(&another) == nullptr) {
    return false;
  }
  const Ellipse& other = dynamic_cast<const Ellipse&>(another);
  return Equal(big_semi_axis_, other.big_semi_axis_) &&
         Equal(small_semi_axis_, other.small_semi_axis_);
}

bool Ellipse::isSimilarTo(const Shape& another) const {
  if (dynamic_cast<const Ellipse*>(&another) == nullptr) {
    return false;
  }
  const Ellipse& other = dynamic_cast<const Ellipse&>(another);
  return Equal(eccentricity(), other.eccentricity());
}

bool Ellipse::containsPoint(const Point& point) const {
  double length_1 = FlatVector(point, focus_1_).Length();
  double length_2 = FlatVector(point, focus_2_).Length();

  return length_1 + length_2 <= 2 * distance_;
}

void Ellipse::rotate(const Point& center, double angle) {
  focus_1_ = focus_1_.rotate(center, angle);
  focus_2_ = focus_2_.rotate(center, angle);
}

void Ellipse::reflex(const Point& center) {
  focus_1_ = focus_1_.reflex(center);
  focus_2_ = focus_2_.reflex(center);
}


void Ellipse::reflex(const Line& axis) {
  focus_1_ = focus_1_.reflex(axis);
  focus_2_ = focus_2_.reflex(axis);
}

void Ellipse::scale(const Point& center, double coefficient) {
  focus_1_ = focus_1_.scale(center, coefficient);
  focus_2_ = focus_2_.scale(center, coefficient);
  distance_ *= coefficient;
  small_semi_axis_ *= coefficient;
  big_semi_axis_ *= coefficient;
}


Circle::Circle(const Point& center, double radius): Ellipse(center, center, 2 * radius) {}

double Circle::radius() const {
  return distance_;
}

Rectangle::Rectangle(const Point& first, const Point& second, double ratio): Polygon(
        std::vector<Point>()) {
  points_.resize(4);
  points_[0] = first;
  points_[2] = second;
  if (ratio > 1) {
    ratio = 1 / ratio;
  }
  points_[1].x = (first.x + ratio * ratio * second.x) / (1 + ratio * ratio) -
                 (second.y - first.y) * ratio / (ratio * ratio + 1);
  points_[1].y = (first.y + ratio * ratio * second.y) / (1 + ratio * ratio) +
                 (second.x - first.x) * ratio / (1 + ratio * ratio);
  points_[3].x = first.x + second.x - points_[1].x;
  points_[3].y = first.y + second.y - points_[1].y;
}

Point Rectangle::center() {
  return Center();
}

std::pair<Line, Line> Rectangle::diagonals() {
  return {Line(points_[0], points_[2]), Line(points_[1], points_[3])};
}

Square::Square(const Point& first, const Point& second): Rectangle(first, second, 1) {}

Circle Square::circumscribedCircle() {
  return Circle(center(), FlatVector(points_[0], points_[2]).Length() / 2);
}


Circle Square::inscribedCircle() {
  return Circle(center(), FlatVector(points_[0], points_[1]).Length() / 2);
}

Triangle::Triangle(const Point& first, const Point& second, const Point& third): Polygon(
        {first, second, third}) {}

Circle Triangle::circumscribedCircle() {
  double length_0 = FlatVector(points_[1], points_[2]).Length();
  double length_1 = FlatVector(points_[0], points_[2]).Length();
  double length_2 = FlatVector(points_[1], points_[0]).Length();
  Point centroid(Center());
  Point orthocenter = this->orthocenter();
  Point center((centroid.x * 3 - orthocenter.x) / 2, (centroid.y * 3 - orthocenter.y) / 2);
  double radius = length_0 * length_1 * length_2 / 4 / area();
  return Circle(center, radius);
}

Circle Triangle::inscribedCircle() {
  double length_0 = FlatVector(points_[1], points_[2]).Length();
  double length_1 = FlatVector(points_[0], points_[2]).Length();
  double length_2 = FlatVector(points_[1], points_[0]).Length();
  Point center((points_[0].x * length_0 + points_[1].x * length_1 + points_[2].x * length_2) /
               (length_0 + length_1 + length_2),
               (points_[0].y * length_0 + points_[1].y * length_1 + points_[2].y * length_2) /
               (length_0 + length_1 + length_2));
  return Circle(center, 2 * area() / perimeter());
}

Point Triangle::centroid() {
  return Center();
}

Point Triangle::orthocenter() {
  Line line_1(points_[0], points_[0].reflex(Line(points_[1], points_[2])));
  Line line_2(points_[1], points_[1].reflex(Line(points_[0], points_[2])));
  return Intersection(line_2, line_1);
}

Line Triangle::EulerLine() {
  return Line(centroid(), orthocenter());
}

Circle Triangle::ninePointsCircle() {
  Circle a = circumscribedCircle();
  a.scale(centroid(), -0.5);
  return a;
}
