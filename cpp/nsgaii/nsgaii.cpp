#include "nsgaii.hpp"
#include <algorithm>
#include <random>
#include <string>
#include <utility>
#include <vector>
#include <cfloat>

namespace NSGAII_0X{

NSGAIIAble::~NSGAIIAble(){}

static std::string log_path = "./NSGAII.log";

static double inf = DBL_MAX;
// 使用dis(gen)即可得到double型均匀分布在[0.0, 1.0]间的随机数
// gain seed for random device
static std::random_device rd;                  
// use rd() to constructe mersenne_twister_engine
static std::mt19937 gen(rd());
// Init random number generator
static std::uniform_real_distribution<double> dis(0.0, 1.0);

struct Individual{
  std::vector<std::vector<char> > genes;// 该个体的基因序列
  std::vector<double> vals;// 优化问题的输入值
  std::vector<double> objs;// 优化问题的目标值
  std::vector<double> constr;// 违反优化问题约束的程度
  double constr_violence;// 总体违反优化问题的约束值
  int rank;// 该个体的被支配等级
  double crowd_distance;// 该个体的拥挤距离

  bool operator<(const Individual& i);
};// struct Individual

// 重载 < 算符, 比较依据为个体的支配等级及拥挤距离
bool Individual::operator<(const Individual& i){
  if (rank < i.rank){
    return true;
  }else if (rank == i.rank && crowd_distance > i.crowd_distance){
    return true;
  }else{
    return false;
  }// if (rank < i.rank)
}// bool Individual::operator <(const Individual& i)

// 用于确认两个个体间的支配关系
//
// 没有违反约束的个体支配违反约束的个体; 
// 都违反约束的情况下通过违法程度来判断, 违反程度小的个体支配违法程度大的个体
// 都未违反约束的情况下通过目标值来判断, 每个目标值都更优的个体处于支配地位
// 返回值  1 表示a支配b
//       -1 表示b支配a
//        0 表示两者间没有支配关系
int CheckDominance(Individual a, Individual b){
  if (a.constr_violence < 0 && b.constr_violence < 0){
    // 两个个体都违反了约束
    if (a.constr_violence > b.constr_violence){
      return 1;
    }else if (a.constr_violence == b.constr_violence){
      return 0;
    }else{
      return -1;
    } // if (a.constr_violence > b.constr_violence)
  }else{
    if (a.constr_violence < 0){
      return -1;
    }else if (b.constr_violence < 0){
      return 1;
    }else{
      // 两个个体都没有违反约束
      int flag_a = 0, flag_b = 0;
      for (size_t i = 0; i < a.objs.size(); i++){
        if (a.objs[i] < b.objs[i]){
          flag_a = 1;
        }
        if (a.objs[i] > b.objs[i]){
          flag_b = 1;
        }
      } // for (int i = 0; i < a.objs.size(); i++)
      if (flag_a == 1 && flag_b == 0){
        return 1;
      }else if (flag_a == 0 && flag_b == 1){
        return -1;
      }else{
        return 0;
      } // if (flag_a == 1 && flag_b == 0)
    } // if (a.constr_violence < 0)
  } // (a.constr_violence < 0 && b.constr_violence < 0)
} // int CheckDominance(Individual a, Individual b)

class NSGAII : NSGAIIAble{
 public:
  explicit NSGAII(Param p);
  ~NSGAII();
  NSGAII(const NSGAII &);
  const NSGAII &operator=(const NSGAII &);

  const Param &GetParam() const{ return this->param; }
  double GetIndividualPropertiesValueByIndex(size_t individual, size_t index) const{ return this->inds_temp[individual].vals[index]; }
  double GetParentPropertiesValue(size_t ind, size_t index) const{ return this->inds[ind].vals[index]; }
  double GetParentObjectValue(size_t ind, size_t index) const{
    return this->inds[ind].objs[index];
  }

  int SetPropertyBoundByIndex(size_t index, double lower_bound, double upper_bound);
  int SetIndividualObjectValueByIndex(size_t individual, size_t index, double value);
  int SetIndividualConstraint(size_t individual, double value);

  int Init();
  // 启动优化
  int EvoluteOnce();

  static NSGAIIAble *GetInstance(Param p) { return new NSGAII(p); }

 private:
  Param param;
  std::vector<Individual> inds;
  std::vector<Individual> inds_temp;
  std::vector<std::pair<double, double> > bound;// 参数区间限制
  size_t scale;

  void non_dominated_sorting();
  void init_inds();
  void decode();
  void crossover();
  void mutate();

  void log(Individual i);
};// class NSGAII

NSGAII::NSGAII(Param p) : param(p){
  // 计算精度值  
  int cnt = p.scale;
  this->scale = 1;
  while (cnt != 0){
    this->scale = this->scale * 10;
    cnt--;
  }// while (cnt != 0)

  for (size_t i = 0; i < p.property_size; ++i){
    bound.push_back(std::make_pair(0.0, 0.0));
  }
}

NSGAII::~NSGAII(){}

// 将临时种群的DNA解码为输入参数
//
// 由于新个体都将在_temp中生成, 故只需对_temp种群解码, 解码公式如下
//
// Positive: if gene[0]==0: value = \sum_{i}^{n} 2^{n-i}*gene[i]
// Negative: if gene[0]==1: value = -(1+\sum_{i}^{n} 2^{n-i}*(!gene[i]))
// value = value / scale
//
void NSGAII::decode(){
  for (size_t ind = 0; ind < this->inds_temp.size(); ++ind){
    this->inds_temp[ind].vals.clear();
    for (size_t i = 0; i < this->inds_temp[ind].genes.size(); ++i){
      double value = 0.0;

      if (this->inds_temp[ind].genes[i][0] == 0){
        // Is positive number
        for (size_t j = 1; j < this->inds_temp[ind].genes[i].size(); ++j){
          value = value * 2 + this->inds_temp[ind].genes[i][j];
        }
      }else{
        // Is negative number
        for (size_t j = 1; j < this->inds_temp[ind].genes[i].size(); j++){
          int g = (this->inds_temp[ind].genes[i][j] == 0)? 1: 0;
          value = value * 2 + g;
        }
        value = -value - 1;
      }// if (this->_temp[ind].genes[i][0] == 0)

      value = value / this->scale;
      // lower than low constrain
      value = (value < this->bound[i].first)? this->bound[i].first: value;
      // higher than high constrain
      value = (value > this->bound[i].second)? this->bound[i].second: value;
      this->inds_temp[ind].vals.push_back(value);
    }// for (size_t i = 0; i < this->inds_temp[ind].genes.size(); ++i)
    
    // 清空目标值
    this->inds_temp[ind].objs.clear();
    for (size_t i = 0; i < this->param.object_size; ++i){
      this->inds_temp[ind].objs.push_back(0.0);
    }
  }// for (size_t ind = 0; ind < this->inds_temp.size(); ++ind)
}// void NSGAII::decode()

// 种群交叉操作
//
// 根据现有的种群产生同等大小的子种群,由于个体已经经过非支配排序,故直接取
// 相邻的个体进行交叉变异
void NSGAII::crossover(){
  for (size_t i = 0; i < this->param.population_size / 2; ++i){
    Individual cind = this->inds[2 * i];
    Individual pind = this->inds[2 * i + 1];
    for (size_t j = 0; j < cind.genes.size(); ++j){
      if (dis(gen) <= this->param.cross_rate){
        size_t point = (size_t)(dis(gen) * (cind.genes[j].size() - 1));
        for (; point < cind.genes[j].size(); ++point){
          std::swap(cind.genes[j][point], pind.genes[j][point]);
        }
      }
    }// for (size_t j = 0; j < cind.genes.size(); ++j)
    this->inds_temp.push_back(cind);
    this->inds_temp.push_back(pind);
  }// for (size_t i = 0; i < this->param.population_size / 2; ++i)
}// void NSGAII::crossover()

// 种群变异操作
// 
// 对新生产的个体进行变异操作
void NSGAII::mutate(){
  for (size_t i = 0; i < this->param.population_size; ++i){
    for (size_t j = 0; j < this->inds[0].genes.size(); ++j){
      if (dis(gen) <= this->param.mutate_rate){
        size_t point = (size_t)(dis(gen) * (this->inds_temp[i].genes[j].size() - 1));
        if (this->inds_temp[i].genes[j][point] == 0){
          this->inds_temp[i].genes[j][point] = 1;
        }else{
          this->inds_temp[i].genes[j][point] = 0;
        }// if (this->inds_temp[i].genes[j][point] == 0)
      }// if (dis(gen) <= this->param.mutate_rate)
    }// for (size_t j = 0; j < this->inds[0].genes.size(); ++j)
  }// for (size_t i = 0; i < this->param.population_size; ++i)
}// void NSGAII::mutate()

void NSGAII::init_inds(){
  // 确定基因链长度, 为了表示负数基因需多一位作为标志位
  std::vector<int> bits;
  double ub, lb, range;
  for (size_t i = 0; i < this->bound.size(); ++i){
    lb = this->bound[i].first;
    ub = this->bound[i].second;
    lb  = (lb  < 0)? -lb : lb;
    ub = (ub < 0)? -ub: ub;
    range = (lb < ub)? ub: lb;
    range = range * scale;
    
    int bit = 2;
    while (range > 1){
      bit++;
      range /= 2;
    }

    bits.push_back(bit);
  }// for (size_t i = 0; i < this->bound.size(); ++i)

  // 初始化inds_temp种群
  for (size_t i = 0; i < this->param.population_size; ++i){
    Individual ind;
    for (size_t j = 0; j < bits.size(); ++j){
      std::vector<char> gene;
      ind.genes.push_back(gene);
      for (int k = 0; k < bits[j]; k++){
        if (dis(gen) > 0.5){
          ind.genes[j].push_back(1);
        }else{
          ind.genes[j].push_back(0);
        }
      } // for (int k = 0; k < bits[j]; k++)
    } // for (int j = 0; j < bits.size(); j++)
    ind.rank = 0;
    ind.crowd_distance = 0.0;
    ind.constr_violence = 0.0;
    this->inds_temp.push_back(ind);
  } // for (int i = 0; i < this->popsize; i++)

  // 解码inds_temp种群
  this->decode();
}// void NSGAII::init_inds()

// 用于对种群进行非支配排序
void NSGAII::non_dominated_sorting(){
  // 储存每个个体支配个体的序号
  std::vector<std::vector<size_t> > dominatedSet;
  // 储存每个个体的被支配数
  std::vector<size_t> dominated;
  for (size_t i = 0; i < this->inds_temp.size(); ++i){
    std::vector<size_t> dSet;
    dominatedSet.push_back(dSet);
    dominated.push_back(0);
  }

  // 计算支配等级
  std::vector<size_t> front;
  for (size_t i = 0; i < this->inds_temp.size(); ++i){
    for (size_t j = i + 1; j < this->inds_temp.size(); ++j){
      int flag = CheckDominance(this->inds_temp[i], this->inds_temp[j]);
      if (flag == 1){
        // inds[i]支配inds[j]
        dominatedSet[i].push_back(j);
        dominated[j]++;
      }else if (flag == -1){
        // inds[j]支配inds[i]
        dominatedSet[j].push_back(i);
        dominated[i]++;
      } 
    }// for (size_t j = i + 1; j < this->inds_temp.size(); ++j)

    if (dominated[i] == 0){
      this->inds_temp[i].rank = 1;
      front.push_back(i);
    } 
  }// for (size_t i = 0; i < this->inds_temp.size(); ++i)

  int rank = 2;
  while(!front.empty()){
    std::vector<size_t> q;
    for (size_t i = 0; i < front.size(); ++i){
      for (size_t j = 0; j < dominatedSet[front[i]].size(); ++j){
        dominated[dominatedSet[front[i]][j]]--;
        if (dominated[dominatedSet[front[i]][j]] == 0){
          this->inds_temp[dominatedSet[front[i]][j]].rank = rank;
          q.push_back(dominatedSet[front[i]][j]);
        }
      }// for (size_t j = 0; j < dominatedSet[front[i]].size(); ++j)
    }// for (size_t i = 0; i < front.size(); ++i)

    rank++;
    front = q;
  }// while(!front.empty())

  // 计算拥挤距离
  std::vector<size_t> index;
  for (size_t i = 0; i < this->inds_temp.size(); ++i){
    index.push_back(i);
  }

  for (size_t i = 0; i < this->inds_temp[0].objs.size(); ++i){
    std::sort(index.begin(), index.end(), 
      [this, i](int a, int b) -> bool {return this->inds_temp[a].objs[i] < this->inds_temp[b].objs[i];});
    this->inds_temp[index[0]].crowd_distance = inf;
    this->inds_temp[index[index.size() - 1]].crowd_distance = inf;
    int min = this->inds_temp[index[0]].objs[i];
    int max = this->inds_temp[index[index.size() - 1]].objs[i];
    int unit = max - min;
    for (size_t j = 1; j < index.size() - 1; ++j){
      if (this->inds_temp[index[j]].crowd_distance != inf){
        this->inds_temp[index[j]].crowd_distance 
        += (this->inds_temp[index[j + 1]].objs[i] - this->inds_temp[index[j - 1]].objs[i]) / unit;
      }
    }// for (size_t j = 1; j < index.size() - 1; ++j)
  }// for (size_t i = 0; i < this->inds_temp[0].objs.size(); ++i)
  
  std::sort(this->inds_temp.begin(), this->inds_temp.end());
}// void NSGAII::non_dominated_sorting()

int NSGAII::SetPropertyBoundByIndex(size_t index, double lower_bound, double upper_bound){
  if (index < 0 || index >= this->param.property_size){
    return 0;
  }
  this->bound[index].first = lower_bound;
  this->bound[index].second = upper_bound;
  return 1;
}

int NSGAII::SetIndividualObjectValueByIndex(size_t individual, size_t index, double value){
  if (individual < 0 || individual >= this->param.population_size){
    return 0;
  }
  if (index < 0 || index >= this->param.object_size){
    return 0;
  }
  this->inds_temp[individual].objs[index] = value;
  return 1;
}

int NSGAII::SetIndividualConstraint(size_t individual, double value){
  if (individual < 0 || individual >= this->param.population_size){
    return 0;
  }
  this->inds_temp[individual].constr.push_back(value);
  return 1;
}

int NSGAII::Init(){
  this->init_inds();

  return 1;
}

int NSGAII::EvoluteOnce(){
  for (size_t i = 0; i < this->inds_temp.size(); ++i){
    this->inds_temp[i].constr_violence = 0;
    for (size_t j = 0; j < this->inds_temp[i].constr.size(); ++j){
      this->inds_temp[i].constr_violence += this->inds_temp[i].constr[j];
    }
    this->inds_temp[i].constr.clear();
  }
  for (size_t i = 0; i < this->inds.size(); ++i){
    Individual ind = this->inds[i];
    this->inds_temp.push_back(ind);
  }
  this->non_dominated_sorting();
  this->inds.clear();
  for (size_t i = 0; i < this->param.population_size; ++i){
    Individual ind = this->inds_temp[i];
    this->inds.push_back(ind);
  }
  this->inds_temp.clear();

  // 生成子代
  this->crossover();
  this->mutate();
  this->decode();

  return 1;
}

void NSGAII::log(Individual i){

}

NSGAIIAble *GetInstance(Param p) { return (NSGAIIAble *)new NSGAII(p); }

/*
#define DLLEXPORT __declspec(dllexport)
extern "C" DLLEXPORT NSGAII *GetInstance(Param p){
  return new NSGAII(p);
}
extern "C" DLLEXPORT NSGAII *GetNSGAII(
  size_t property_size, 
  size_t object_size, 
  size_t population_size=16, 
  double cross_rate=0.6, 
  double mutate_rate=0.2,
  int scale=3){
  Param p(population_size, object_size, population_size, cross_rate, mutate_rate, scale);
  return new NSGAII(p);
}
extern "C" DLLEXPORT double GetIndividualPropertiesValueByIndex(NSGAIIAble *ptr, size_t individual, size_t index){
  return ptr->GetIndividualPropertiesValueByIndex(individual, index);
}
extern "C" DLLEXPORT double GetParentPropertiesValue(NSGAIIAble *ptr, size_t ind, size_t index){
  return ptr->GetParentPropertiesValue(ind, index);
}
extern "C" DLLEXPORT double GetParentObjectValue(NSGAIIAble *ptr, size_t ind, size_t index){
  return ptr->GetParentObjectValue(ind, index);
}
extern "C" DLLEXPORT int SetPropertyBoundByIndex(NSGAIIAble *ptr, size_t index, double lower_bound, double upper_bound){
  return ptr->SetPropertyBoundByIndex(index, lower_bound, upper_bound);
}
extern "C" DLLEXPORT int SetIndividualObjectValueByIndex(NSGAIIAble *ptr, size_t individual, size_t index, double value){
  return ptr->SetIndividualObjectValueByIndex(individual, index, value);
}
extern "C" DLLEXPORT int SetIndividialConstraint(NSGAIIAble *ptr, size_t individual, double value){
  return ptr->SetIndividualConstraint(individual, value);
}
extern "C" DLLEXPORT int Init(NSGAIIAble *ptr){
  return ptr->Init();
}
extern "C" DLLEXPORT int EvoluteOnce(NSGAIIAble *ptr){
  return ptr->EvoluteOnce();
}
#undef DLLEXPORT
*/
}// NSGAII_0X
