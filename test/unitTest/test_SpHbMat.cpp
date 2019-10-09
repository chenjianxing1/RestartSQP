#include <sqphot/Matrix.hpp>

template <typename T>
bool TEST_EQUAL(T a, T b){
    return a == b;
}

bool TEST_EQUAL_INT_ARRAY(int* a, int* b, int length, const char* name = NULL){
    for(int i = 0; i < length; i++){
        if(a[i] != b[i]){
            printf("TEST_EQUAL_INT_ARRAY for %s failed at index %d\n", name, i);
            return false;
        }
    }
    return true;
    
}

bool TEST_EQUAL_DOUBLE_ARRAY(double* a, double* b, int length, const char* name = NULL){
    for(int i = 0; i < length; i++){
        if(a[i] != b[i]){
            printf("TEST_EQUAL_DOUBLE_ARRAY for %s failed at index %d\n", name, i);
            return false;
        }
    }
    return true;
}
using namespace SQPhotstart;
int main(int argc, char* argv[]){


double data[6]= {1,2,3,4,5,6};
int rowIndex_csr[3] = {0,3,6};
int colIndex_csr[6] = {0,1,2,0,1, 2};
double matval_csr[6] = {1,2,3,4,5,6}; 

int colIndex_csc[4] = {0,2,4,6};
int rowIndex_csc[6] = {0,1,0,1,0,1};
double matval_csc[6] ={1,4,2,5,3,6};

//shared_ptr<SpHbMat> m_csr_col_oriented = make_shared<SpHbMat>(data,2,3,false,false);
//shared_ptr<SpHbMat> m_csc_col_oriented = make_shared<SpHbMat>(data,2,3,true,true);

shared_ptr<SpHbMat> m_csr_row_oriented = make_shared<SpHbMat>(data,2,3,true,true);
shared_ptr<SpHbMat> m_csc_row_oriented = make_shared<SpHbMat>(data,2,3,true,false);
TEST_EQUAL_INT_ARRAY(m_csc_row_oriented->RowIndex(),rowIndex_csc,6,"csc_row");
TEST_EQUAL_INT_ARRAY(m_csc_row_oriented->ColIndex(),colIndex_csc,4,"csc_col");
TEST_EQUAL_DOUBLE_ARRAY(m_csc_row_oriented->MatVal(),matval_csc,6,"csc_mat");
TEST_EQUAL_INT_ARRAY(m_csr_row_oriented->RowIndex(),rowIndex_csr,3,"csr_row");
TEST_EQUAL_INT_ARRAY(m_csr_row_oriented->ColIndex(),colIndex_csr,6,"csr_col");
TEST_EQUAL_DOUBLE_ARRAY(m_csr_row_oriented->MatVal(),matval_csr,6,"csr_mat");
printf("---------------------------------------------------------\n");
printf("    Converting dense matrix to Sparse matrix succeeded!   \n");
printf("---------------------------------------------------------\n");

}

