/**
 * @file   matcher.h
 * @author Ruo Li <rli@aztec>
 * @date   Sun Apr  4 17:53:15 2010
 * 
 * @brief  将 (-x, 0) 和 (x, 0) 进行匹配的匹配器。
 * 
 * 
 */

struct my_matcher {
  int bmark;
  /** 
   * 对点 p0 和 p1 进行坐标匹配。
   */
  int value(const Point<2>& p0,
            const Point<2>& p1,
            double tol) const {
    double d0 = fabs(p0[0] - p1[0]);
    double d1 = fabs(p0[1] - p1[1]);
    if (d0 + d1 < tol) return 0;
    d0 = fabs(p0[0] + p1[0]);
    if (d0 + d1 < tol) return 1;
    return -1;
  }
  /** 
   * 读入周期描述的配置文件。
   */
  void readConfigFile(const std::string& file) {
    filtering_istream is;
    OpenFilteredStream(file, is);
    is >> bmark;
  }
};

/**
 * end of file
 * 
 */
