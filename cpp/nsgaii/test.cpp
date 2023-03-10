#include "nsgaii.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
using std::cin;
using std::cout;
using std::string;
using std::vector;

int main() {
    int cnt = 2;
    int ind_cnt = 48, obj_cnt = 4;
    const double constr[] = {160, -30, 0, 0};
    NSGAII_0X::Param p = NSGAII_0X::Param(1, obj_cnt, ind_cnt);
    NSGAII_0X::NSGAIIAble *ptr = NSGAII_0X::GetInstance(p); ptr->Init();
    for (int pro = 0; pro < 1; ++pro) { ptr->SetPropertyBoundByIndex(pro, -1, 1); }

    vector<vector<double>> res(ind_cnt, vector<double>(obj_cnt, 0));

    while (cnt) {
        string tmp = "";
        while (tmp[0] != '/') { std::getline(cin, tmp); }
        for (int ind = 0; ind < ind_cnt; ++ind) {
            for (int obj = 0; obj < obj_cnt; ++obj) {
                double v; cin >> v;
                ptr->SetIndividualObjectValueByIndex(ind, obj, -1*v);
                ptr->SetIndividualConstraint(ind, std::min(0.0, v-constr[obj]));
            }
        }
        ptr->EvoluteOnce();

        for (int ind = 0; ind < ind_cnt; ++ind) {
            for (int obj = 0; obj < obj_cnt; ++obj) {
                res[ind][obj] = -1 * ptr->GetParentObjectValue(ind, obj);
            }
        }

        std::sort(res.begin(), res.end(), [](vector<double> a, vector<double> b) -> bool {
            bool flag = true;
            for (size_t i = 0; i < a.size(); ++i) {
                if (a[i] <= b[i]) { flag = false; break; }
            }
            return flag;
        });
        
        for (auto r : res) {
            for (auto o : r) { cout << o << " "; }
            cout << "\n";
        }
        cout << "\n";

        --cnt;
    }

    return 0;
}