#include "nsgaii.hpp"
#include <iostream>
#include <string>
using std::cin;
using std::cout;
using std::string;

int main() {
    int cnt = 2;
    int ind_cnt = 48, obj_cnt = 4;
    NSGAII_0X::Param p = NSGAII_0X::Param(1, obj_cnt, ind_cnt);
    NSGAII_0X::NSGAIIAble *ptr = NSGAII_0X::GetInstance(p); ptr->Init();
    for (int pro = 0; pro < 1; ++pro) { ptr->SetPropertyBoundByIndex(pro, -1, 1); }

    while (cnt) {
        string tmp = "";
        while (tmp[0] != '/') { std::getline(cin, tmp); }
        for (int ind = 0; ind < ind_cnt; ++ind) {
            for (int obj = 0; obj < obj_cnt; ++obj) {
                double v; cin >> v;
                ptr->SetIndividualObjectValueByIndex(ind, obj, -1*v);
            }
        }
        ptr->EvoluteOnce();

        for (int ind = 0; ind < ind_cnt; ++ind) {
            cout << ind << ": ";
            for (int obj = 0; obj < obj_cnt; ++obj) {
                cout << -1 * ptr->GetParentObjectValue(ind, obj) << " "; 
            }
            cout << "\n";
        }
        cout << "\n";
        --cnt;
    }

    return 0;
}