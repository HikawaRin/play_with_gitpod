#ifndef NSGAII_0X_SRC_NSGAII_HPP
#define NSGAII_0X_SRC_NSGAII_HPP
#include <climits>
using std::size_t;

namespace NSGAII_0X
{
// NSGAII算法接口, 暂时仅提供Windwos动态链接库
//
// 在加载DLL后通过 NSGAIIAble *nsgaii_ptr = NSGAII_0X::GetInstance(param) 来获取一个算法实体
// 在获得算法实体后调用 SetPropertyBoundByIndex 设置每个自变量的范围
// 
// 第一次调用 nsgaii_ptr->Evolution() 算法会自动完成初始化
// 通过 nsgaii_ptr->Evolution() 可以进行一次进化, 在每次进化前确保每个个体的目标值已被计算且写入对应位置

struct Param{
  size_t population_size;// 种群数量
  double cross_rate;// 交叉概率
  double mutate_rate;// 变异概率
  size_t property_size;// 自变量数量
  size_t object_size;// 目标数量
  int scale;// 参数精度

  Param(size_t property_size, 
  size_t object_size, 
  size_t population_size=16, 
  double cross_rate=0.6, 
  double mutate_rate=0.2,
  int scale=3) 
  : population_size(population_size), 
  cross_rate(cross_rate), 
  mutate_rate(mutate_rate), 
  property_size(property_size), 
  object_size(object_size),
  scale(scale){}
};

class NSGAIIAble{
 public:
  virtual ~NSGAIIAble()=0;
  virtual const Param &GetParam() const=0;
  virtual double GetIndividualPropertiesValueByIndex(size_t individual, size_t index) const=0;
  virtual double GetParentPropertiesValue(size_t ind, size_t index) const=0;
  virtual double GetParentObjectValue(size_t ind, size_t index) const=0;

  virtual int SetPropertyBoundByIndex(size_t index, double lower_bound, double upper_bound)=0;
  virtual int SetIndividualObjectValueByIndex(size_t individual, size_t index, double value)=0;
  virtual int SetIndividualConstraint(size_t individual, double value)=0;

  virtual int Init()=0;
  // 进化一次
  // 在进行进化前需要写入所有个体的目标值
  virtual int EvoluteOnce()=0;
  // TODO 
  // 通过该函数可以判断此时种群是否收敛
  // 
  // virtual bool IsConverge()=0;

  static NSGAIIAble *GetInstance(Param p);
};// class NSGAII

} // namespace NSGAII_0X
#endif// NSGAII_0X_SRC_NSGAII_HPP
