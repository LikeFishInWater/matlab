#include "mex.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

//[qAry, TABLE_MULTIPLY, TABLE_ADD, TABLE_INVERSE, maxVarDegree, maxChkDegree, VarDegree, ChkDegree,H, Row_Link_Col, Col_Link_Row,Q, R, L_Post, LDR_Vector, L_SIGMA, L_RHO,CodeLen, ChkLen, MsgLen,Rate] = Initial_CPP('PEG.64.128.gf.8.txt.','Table.GF.8.txt');

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    if(nrhs!=2) mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "2 input parameters required.");
    if(nlhs!=21) mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "21 output parameters required.");
    
    string CodeFileName = mxArrayToString(prhs[0]);
    string TableFileName = mxArrayToString(prhs[1]);
    
    // 输出的变量
    unsigned int qAry, **TABLE_MULTIPLY, **TABLE_ADD, *TABLE_INVERSE;   // 0~3
    unsigned int maxVarDegree, maxChkDegree, *VarDegree, *ChkDegree;     // 4~7
    unsigned int **H;   // 8
    unsigned long **Row_Link_Col, **Col_Link_Row;     // 9~10
    double ***Q, ***R, ** L_Post, * LDR_Vector, * L_SIGMA, * L_RHO;   // 11~16
    unsigned long CodeLen, ChkLen, MsgLen;    // 17~19
    double Rate;   // 20
    
    ifstream fin(CodeFileName);
    if (!fin.is_open()) mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Can't open CodeFile");
    
    int tempqAry;
    fin >> CodeLen >> ChkLen >> tempqAry;
    qAry = unsigned int(tempqAry);
    MsgLen = CodeLen - ChkLen;
    Rate = double(MsgLen) / double(CodeLen);//R=K/N
    fin >> maxVarDegree >> maxChkDegree;
    
    VarDegree = new unsigned int[CodeLen];
    ChkDegree = new unsigned int[ChkLen];
    Col_Link_Row = new unsigned long*[CodeLen];
    Row_Link_Col = new unsigned long*[ChkLen];
    for (unsigned long i = 0; i < CodeLen; i++) { fin >> VarDegree[i]; Col_Link_Row[i] = new unsigned long[maxVarDegree]; }
    for (unsigned long j = 0; j < ChkLen; j++) { fin >> ChkDegree[j]; Row_Link_Col[j] = new unsigned long[maxChkDegree]; }
    
    H = new unsigned int*[ChkLen]();
    for (unsigned long j = 0; j < ChkLen; j++) H[j] = new unsigned int[CodeLen]();
    
    unsigned long tempIndex; unsigned int tempQARY;
    for (unsigned long i = 0; i < CodeLen; i++)
    {
        for (unsigned int di = 0; di < VarDegree[i]; di++)
        {
            fin >> tempIndex >> tempQARY;
            if (tempIndex > 0)
            {
                Col_Link_Row[i][di] = tempIndex - 1;
                H[tempIndex - 1][i] = unsigned int(tempQARY);
            }
        }
    }
    
    for (unsigned long j = 0; j < ChkLen; j++)
    {
        for (unsigned int dj = 0; dj < ChkDegree[j]; dj++)
        {
            fin >> tempIndex >> tempQARY;
            if (tempIndex > 0)
            {
                Row_Link_Col[j][dj] = tempIndex - 1;
                H[j][tempIndex - 1] = tempQARY;
            }
        }
    }
    
    fin.close();
    
    int tempValue;
    ifstream fin_table(TableFileName);
    if (!fin_table.is_open()) mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Can't open TableFile");
    
    string rub;
    fin_table >> rub >> rub;
    
    TABLE_MULTIPLY = new unsigned int* [qAry];
    for (unsigned int i = 0; i < qAry; i++)
    {
        TABLE_MULTIPLY[i] = new unsigned int[qAry];
        for(unsigned int j = 0; j < qAry; j++)
        {
            fin_table >> tempValue;
            TABLE_MULTIPLY[i][j] = unsigned int(tempValue);
        }
    }
    
    fin_table >> rub >> rub;
    TABLE_INVERSE = new unsigned int [qAry];
    TABLE_INVERSE[0] = 0;
    for (unsigned int k = 1; k < qAry; k++)
    {
        fin_table >> tempValue;
        TABLE_INVERSE[k] = unsigned int(tempValue);
    }
    
    TABLE_ADD = new unsigned int* [qAry];
    for (unsigned int i = 0; i < qAry; i++)
    {
        TABLE_ADD[i] = new unsigned int[qAry];
        for (unsigned int j = 0; j < qAry; j++){
            TABLE_ADD[i][j] = i ^ j;
        }
    }
    
    Q = new double **[ChkLen]();
    R = new double **[ChkLen]();
    for(unsigned long j = 0; j < ChkLen; j++)
    {
        Q[j] = new double * [CodeLen]();
        R[j] = new double * [CodeLen]();
        for(unsigned long i = 0; i < CodeLen; i++)
        {
            Q[j][i] = new double[qAry - 1];
            R[j][i] = new double[qAry - 1];
            for(int k=0;k<qAry-1;k++){
                Q[j][i][k]= 0;
                R[j][i][k]=0;
            }
        }
    }
    
    L_Post = new double*[CodeLen]();
    for (unsigned long i = 0; i < CodeLen; i++){
        L_Post[i] = new double[qAry - 1]();
        for (int j=0;j<qAry-1;j++)
            L_Post[i][j] = 0;
    }
    
    LDR_Vector = new double[qAry - 1]; for (int j=0;j<qAry-1;j++) LDR_Vector[j] = 0;
    L_SIGMA = new double[qAry - 1]; for (int j=0;j<qAry-1;j++) L_SIGMA[j] = 0;
    L_RHO = new double[qAry - 1]; for (int j=0;j<qAry-1;j++) L_RHO[j] = 0;
    
    // 变量复制到输出
    int dims_2[] = {1, 1};
    plhs[17] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS, mxREAL); memcpy(mxGetPr(plhs[17]), &CodeLen, sizeof(unsigned long));
    plhs[18] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS, mxREAL); memcpy(mxGetPr(plhs[18]), &ChkLen, sizeof(unsigned long));
    plhs[0] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS, mxREAL); memcpy(mxGetPr(plhs[0]), &qAry, sizeof(unsigned int));
    plhs[19] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS, mxREAL); memcpy(mxGetPr(plhs[19]), &MsgLen, sizeof(unsigned long));
    plhs[20] = mxCreateDoubleMatrix(1,1, mxREAL); memcpy(mxGetPr(plhs[20]), &Rate, sizeof(double));
    plhs[4] = mxCreateNumericArray(2,dims_2, mxUINT32_CLASS, mxREAL); memcpy(mxGetPr(plhs[4]), &maxVarDegree, sizeof(unsigned int));
    plhs[5] = mxCreateNumericArray(2,dims_2, mxUINT32_CLASS, mxREAL); memcpy(mxGetPr(plhs[5]), &maxChkDegree, sizeof(unsigned int));
    
    dims_2[0] = 1; dims_2[1] = CodeLen;
    plhs[6] = mxCreateNumericArray(2,dims_2, mxUINT32_CLASS, mxREAL); memcpy(mxGetPr(plhs[6]), VarDegree, CodeLen*sizeof(unsigned int));
    dims_2[0] = 1; dims_2[1] = ChkLen;
    plhs[7] = mxCreateNumericArray(2,dims_2, mxUINT32_CLASS, mxREAL); memcpy(mxGetPr(plhs[7]), ChkDegree, ChkLen*sizeof(unsigned int));
    
    dims_2[0] = CodeLen; dims_2[1] = maxVarDegree;
    plhs[10] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS ,mxREAL);
    unsigned long *ptr_ulong = (unsigned long *)mxGetData(plhs[10]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int di = 0; di < dims_2[1]; di++) memcpy(ptr_ulong+i+di*dims_2[0], Col_Link_Row[i]+di, sizeof(unsigned long));
    }
    
    dims_2[0] = ChkLen; dims_2[1] = maxChkDegree;
    plhs[9] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS ,mxREAL);
    ptr_ulong = (unsigned long *)mxGetData(plhs[9]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int di = 0; di < dims_2[1]; di++) memcpy(ptr_ulong+i+di*dims_2[0], Row_Link_Col[i]+di, sizeof(unsigned long));
    }
    
    dims_2[0] = ChkLen; dims_2[1] = CodeLen;
    plhs[8] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS, mxREAL); 
    unsigned int *ptr_uint = (unsigned int *)mxGetData(plhs[8]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int di = 0; di < dims_2[1]; di++) memcpy(ptr_uint+i+di*dims_2[0], H[i]+di, sizeof(unsigned int));
    }
    
    dims_2[0] = qAry; dims_2[1] = qAry;
    plhs[1] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS, mxREAL); 
    ptr_uint = (unsigned int *)mxGetData(plhs[1]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int di = 0; di < dims_2[1]; di++) memcpy(ptr_uint+i+di*dims_2[0], TABLE_MULTIPLY[i]+di, sizeof(unsigned int));
    }
    
    dims_2[0] = 1; dims_2[1] = qAry;
    plhs[3] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetPr(plhs[3]), TABLE_INVERSE, qAry*sizeof(unsigned int));
    
    dims_2[0] = qAry; dims_2[1] = qAry;
    plhs[2] = mxCreateNumericArray(2, dims_2, mxUINT32_CLASS, mxREAL); 
    ptr_uint = (unsigned int *)mxGetData(plhs[2]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int di = 0; di < dims_2[1]; di++) memcpy(ptr_uint+i+di*dims_2[0], TABLE_ADD[i]+di, sizeof(unsigned int));
    }
    
    int dims_3[] = {ChkLen, CodeLen, qAry-1};
    plhs[11] = mxCreateNumericArray(3, dims_3, mxDOUBLE_CLASS, mxREAL);
    double *ptr_double = (double *)mxGetData(plhs[11]);
    for (int i = 0; i < dims_3[0]; i++)
    {
        for (int j = 0; j < dims_3[1]; j++){
            for (int k = 0; k < dims_3[2]; k++)
                memcpy(ptr_double+(k*dims_3[1]+j)*dims_3[0]+i, Q[i][j]+k, sizeof(double));
        }
    }
    
    plhs[12] = mxCreateNumericArray(3, dims_3, mxDOUBLE_CLASS, mxREAL); 
    ptr_double = (double *)mxGetData(plhs[12]);
    for (int i = 0; i < dims_3[0]; i++)
    {
        for (int j = 0; j < dims_3[1]; j++){
            for (int k = 0; k < dims_3[2]; k++)
                memcpy(ptr_double+(k*dims_3[1]+j)*dims_3[0]+i, R[i][j]+k, sizeof(double));
        }
    }
    
    dims_2[0] = CodeLen; dims_2[1] = qAry-1;
    plhs[13] = mxCreateNumericArray(2, dims_2, mxDOUBLE_CLASS, mxREAL); 
    ptr_double = (double *)mxGetData(plhs[13]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int di = 0; di < dims_2[1]; di++)
            memcpy(ptr_double+i+di*dims_2[0], L_Post[i]+di, sizeof(double));
    }
    
    plhs[14] = mxCreateDoubleMatrix(1, qAry - 1, mxREAL);
    ptr_double = (double *)mxGetData(plhs[14]);
    for (int i = 0; i < qAry - 1; i++) memcpy(ptr_double+i, LDR_Vector+i, sizeof(double));

    plhs[15] = mxCreateDoubleMatrix(1, qAry - 1, mxREAL);
    ptr_double = (double *)mxGetData(plhs[15]);
    for (int i = 0; i < qAry - 1; i++) memcpy(ptr_double+i, L_SIGMA+i, sizeof(double));
   
    plhs[16] = mxCreateDoubleMatrix(1, qAry - 1, mxREAL); 
    ptr_double = (double *)mxGetData(plhs[16]);
    for (int i = 0; i < qAry - 1; i++) memcpy(ptr_double+i, L_RHO+i, sizeof(double));
    
    // 释放不需要的变量的空间
    delete[] VarDegree;
    delete[] ChkDegree;
    for (unsigned long i = 0; i < CodeLen; i++) delete[] Col_Link_Row[i]; delete[] Col_Link_Row;
    for (unsigned long i = 0; i < ChkLen; i++) delete[] Row_Link_Col[i]; delete[] Row_Link_Col;
    
    for (unsigned long j = 0; j < ChkLen; j++) delete[] H[j]; delete[] H;
    
    for (unsigned int i = 0; i < qAry; i++) delete[] TABLE_MULTIPLY[i]; delete[] TABLE_MULTIPLY;
    
    delete[] TABLE_INVERSE;
    
    for (unsigned int i = 0; i < qAry; i++) delete[] TABLE_ADD[i]; delete[] TABLE_ADD;
    
    for(unsigned long j = 0; j < ChkLen; j++)
    {
        for(unsigned long i = 0; i < CodeLen; i++)
        {
            delete[] Q[j][i]; delete[] R[j][i];
        }
        delete[] Q[j]; delete[] R[j];
    }
    delete[] Q; delete[] R;
    
    for (unsigned long i = 0; i < CodeLen; i++) delete[] L_Post[i]; delete[] L_Post;
    delete[] LDR_Vector;
    delete[] L_SIGMA;
    delete[] L_RHO;
}