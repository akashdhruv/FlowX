#include <iostream>
#include<stdio.h>

using namespace std;

class unit_interface {
    public:
        virtual void unit_initialize(int *nx, int *ny) { cout <<"Constructing Interface"<<endl;};
};   


class unit_main : public unit_interface{
    public:
        void unit_initialize(int *nx, int *ny) {
            unit_nx = nx;
            unit_ny = ny;

            unit_ptr = unit_ny;
            cout <<"Constructing Main "<<*unit_nx+*unit_ny<<" "<<unit_ptr<<endl;
        }

        int *unit_ny;

    private:
        int *unit_nx;
        int *unit_ptr;

        
};

int main(){

    int nx,ny;
    int *ptr;

    nx = 20;
    ny = 40;

    ptr = &ny;

    unit_main unit;

    unit.unit_initialize(&nx, &ny);

    cout<<"Accessing object data from main "<<*unit.unit_ny<<" "<<ptr<<endl;

    return 0;
}; 
