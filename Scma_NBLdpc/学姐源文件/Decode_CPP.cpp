#include "mex.h"
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

bool Decode(double ** LDR, unsigned int *DecodingSeq,    // 0~1
        unsigned long CodeLen, unsigned long ChkLen,      // 2~3
        unsigned int *VarDegree, unsigned int *ChkDegree,   // 4~5
        unsigned long **Col_Link_Row, unsigned long **Row_Link_Col,   // 6~7
        unsigned int qAry, unsigned int **TABLE_MULTIPLY, unsigned int **TABLE_ADD, unsigned int *TABLE_INVERSE,    // 8~11
        unsigned int **H,    // 12
        double ***Q, double  ***R,  // 13~14
        double  **L_Post, double *LDR_Vector, double *L_SIGMA, double *L_RHO,   // 15~18
        unsigned int maxIteration);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=19) mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "19 input parameters required.");
    if(nlhs!=2) mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "2 output parameters required.");
    
    // 输入的变量
    double ** LDR;   // 0
    unsigned long CodeLen, ChkLen;   // 1~2
    unsigned int *VarDegree, *ChkDegree;   // 3~4
    unsigned long **Col_Link_Row, **Row_Link_Col;   // 5~6
    unsigned int qAry, **TABLE_MULTIPLY, **TABLE_ADD, *TABLE_INVERSE;    // 7~10
    unsigned int **H;    // 11
    double ***Q, ***R;  // 12~13
    double  **L_Post, *LDR_Vector, *L_SIGMA, *L_RHO;   // 14~17
    unsigned int maxIteration; //18
    
    // 输出的变量
    unsigned int *DecodingSeq;    // 0
    double *Decode_output; //1
    
    size_t dims_2[] = {0,0};
    int * dims_3;
    double *ptr_double;
    unsigned long * ptr_ulong;
    unsigned int * ptr_uint;
    
    // 输入变量分配空间，并从matlab的输入参数复制
    // prhs[0]
    dims_2[0] = mxGetM(prhs[0]); dims_2[1] = mxGetN(prhs[0]);
    LDR = new double*[dims_2[0]];
    for (int i=0;i<dims_2[0];i++) LDR[i] = new double[dims_2[1]];
    ptr_double = (double *)mxGetData(prhs[0]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int j = 0; j < dims_2[1]; j++)
            memcpy(LDR[i]+j, ptr_double+i+j*dims_2[0], sizeof(double));
    }
    // prhs[1]
    ptr_ulong = (unsigned long *)mxGetData(prhs[1]);
    memcpy(&CodeLen, ptr_ulong, sizeof(unsigned long));
    // prhs[2]
    ptr_ulong = (unsigned long *)mxGetData(prhs[2]);
    memcpy(&ChkLen, ptr_ulong, sizeof(unsigned long));
    // prhs[3]
    dims_2[0] = mxGetM(prhs[3]); dims_2[1] = mxGetN(prhs[3]);
    VarDegree = new unsigned int[dims_2[0]*dims_2[1]];
    ptr_uint = (unsigned int *)mxGetData(prhs[3]);
    memcpy(VarDegree, ptr_uint, dims_2[0]*dims_2[1]*sizeof(unsigned int));
    // prhs[4]
    dims_2[0] = mxGetM(prhs[4]); dims_2[1] = mxGetN(prhs[4]);
    ChkDegree = new unsigned int[dims_2[0]*dims_2[1]];
    ptr_uint = (unsigned int *)mxGetData(prhs[4]);
    memcpy(ChkDegree, ptr_uint, dims_2[0]*dims_2[1]*sizeof(unsigned int));
    // prhs[5]
    dims_2[0] = mxGetM(prhs[5]); dims_2[1] = mxGetN(prhs[5]);
    Col_Link_Row = new unsigned long*[dims_2[0]];
    for (int i=0;i<dims_2[0];i++) Col_Link_Row[i] = new unsigned long[dims_2[1]];
    ptr_ulong = (unsigned long *)mxGetData(prhs[5]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int j = 0; j < dims_2[1]; j++)
            memcpy(Col_Link_Row[i]+j, ptr_ulong+i+j*dims_2[0], sizeof(unsigned long));
    }
    // prhs[6]
    dims_2[0] = mxGetM(prhs[6]); dims_2[1] = mxGetN(prhs[6]);
    Row_Link_Col = new unsigned long*[dims_2[0]];
    for (int i=0;i<dims_2[0];i++) Row_Link_Col[i] = new unsigned long[dims_2[1]];
    ptr_ulong = (unsigned long *)mxGetData(prhs[6]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int j = 0; j < dims_2[1]; j++)
            memcpy(Row_Link_Col[i]+j, ptr_ulong+i+j*dims_2[0], sizeof(unsigned long));
    }
    // prhs[7]
    ptr_uint = (unsigned int *)mxGetData(prhs[7]);
    memcpy(&qAry, ptr_uint, sizeof(unsigned int));
    // prhs[8]
    dims_2[0] = mxGetM(prhs[8]); dims_2[1] = mxGetN(prhs[8]);
    TABLE_MULTIPLY = new unsigned int*[dims_2[0]];
    for (int i=0;i<dims_2[0];i++) TABLE_MULTIPLY[i] = new unsigned int[dims_2[1]];
    ptr_uint = (unsigned int *)mxGetData(prhs[8]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int j = 0; j < dims_2[1]; j++)
            memcpy(TABLE_MULTIPLY[i]+j, ptr_uint+i+j*dims_2[0], sizeof(unsigned int));
    }
    // prhs[9]
    dims_2[0] = mxGetM(prhs[9]); dims_2[1] = mxGetN(prhs[9]);
    TABLE_ADD = new unsigned int*[dims_2[0]];
    for (int i=0;i<dims_2[0];i++) TABLE_ADD[i] = new unsigned int[dims_2[1]];
    ptr_uint = (unsigned int *)mxGetData(prhs[9]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int j = 0; j < dims_2[1]; j++)
            memcpy(TABLE_ADD[i]+j, ptr_uint+i+j*dims_2[0], sizeof(unsigned int));
    }
    // prhs[10]
    dims_2[0] = mxGetM(prhs[10]); dims_2[1] = mxGetN(prhs[10]);
    TABLE_INVERSE = new unsigned int[dims_2[0]*dims_2[1]];
    ptr_uint = (unsigned int *)mxGetData(prhs[10]);
    memcpy(TABLE_INVERSE, ptr_uint, dims_2[0]*dims_2[1]*sizeof(unsigned int));
    // prhs[11]
    dims_2[0] = mxGetM(prhs[11]); dims_2[1] = mxGetN(prhs[11]);
    H = new unsigned int*[dims_2[0]];
    for (int i=0;i<dims_2[0];i++) H[i] = new unsigned int[dims_2[1]];
    ptr_uint = (unsigned int *)mxGetData(prhs[11]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int j = 0; j < dims_2[1]; j++)
            memcpy(H[i]+j, ptr_uint+i+j*dims_2[0], sizeof(unsigned int));
    }
    
    // prhs[12]
    dims_3 = (int *)mxGetDimensions(prhs[12]);
    Q = new double**[dims_3[0]];
    for (int i=0;i<dims_3[0];i++){
        Q[i] = new double*[dims_3[1]];
        for (int j=0;j<dims_3[1];j++)
            Q[i][j] = new double[dims_3[2]];
    }
    ptr_double = (double *)mxGetData(prhs[12]);
    for (int i = 0; i < dims_3[0]; i++)
    {
        for (int j = 0; j < dims_3[1]; j++){
            for (int k = 0; k < dims_3[2]; k++)
                //memcpy(Q[i][j]+k, ptr_double+i+(k*dims_3[0]+j)*dims_3[1], sizeof(double));
                memcpy(Q[i][j]+k, ptr_double+(k*dims_3[1]+j)*dims_3[0]+i, sizeof(double));
        }
    }
     
    // prhs[13]
    dims_3 = (int *)mxGetDimensions(prhs[13]);
    R = new double**[dims_3[0]];
    for (int i=0;i<dims_3[0];i++){
        R[i] = new double*[dims_3[1]];
        for (int j=0;j<dims_3[1];j++)
            R[i][j] = new double[dims_3[2]];
    }
    
    ptr_double = (double *)mxGetData(prhs[13]);
    for (int i = 0; i < dims_3[0]; i++)
    {
        for (int j = 0; j < dims_3[1]; j++){
            for (int k = 0; k < dims_3[2]; k++)
                memcpy(R[i][j]+k, ptr_double+(k*dims_3[1]+j)*dims_3[0]+i, sizeof(double));
        }
    }
    
    // prhs[14]
    dims_2[0] = mxGetM(prhs[14]); dims_2[1] = mxGetN(prhs[14]);
    L_Post = new double*[dims_2[0]];
    for (int i=0;i<dims_2[0];i++) L_Post[i] = new double[dims_2[1]];
    ptr_double = (double *)mxGetData(prhs[14]);
    for (int i = 0; i < dims_2[0]; i++)
    {
        for (int j = 0; j < dims_2[1]; j++)
            memcpy(L_Post[i]+j, ptr_double+i+j*dims_2[0], sizeof(double));
    }
    
    // prhs[15]
    dims_2[0] = mxGetM(prhs[15]); dims_2[1] = mxGetN(prhs[15]);
    LDR_Vector = new double[dims_2[0]*dims_2[1]];
    ptr_double = (double *)mxGetData(prhs[15]);
    memcpy(LDR_Vector, ptr_double, dims_2[0]*dims_2[1]*sizeof(double));
    
    // prhs[16]
    dims_2[0] = mxGetM(prhs[16]); dims_2[1] = mxGetN(prhs[16]);
    L_SIGMA = new double[dims_2[0]*dims_2[1]];
    ptr_double = (double *)mxGetData(prhs[16]);
    memcpy(L_SIGMA, ptr_double, dims_2[0]*dims_2[1]*sizeof(double));
    
    // prhs[17]
    dims_2[0] = mxGetM(prhs[17]); dims_2[1] = mxGetN(prhs[17]);
    L_RHO = new double[dims_2[0]*dims_2[1]];
    ptr_double = (double *)mxGetData(prhs[17]);
    memcpy(L_RHO, ptr_double, dims_2[0]*dims_2[1]*sizeof(double));

    // prhs[18]
   // ptr_uint = (unsigned int *)mxGetData(prhs[18]);
    //memcpy(&maxIteration, ptr_uint, sizeof(unsigned int));
    ptr_double=mxGetPr(prhs[18]);
	maxIteration=(unsigned int) *ptr_double;
    // 输出变量分配空间
    DecodingSeq = new unsigned int[CodeLen];
    Decode_output = new double[CodeLen];

    // 调用的Cpp解码函数
    Decode(LDR, DecodingSeq, CodeLen, ChkLen, VarDegree, ChkDegree, Col_Link_Row, Row_Link_Col, qAry, TABLE_MULTIPLY, TABLE_ADD, TABLE_INVERSE,
            H, Q, R, L_Post, LDR_Vector, L_SIGMA, L_RHO, maxIteration);
    
    // 释放输入参数分配的空间
    dims_2[0] = mxGetM(prhs[0]);
    for (int i=0;i<dims_2[0];i++) delete[] LDR[i]; delete[] LDR;    // 0
    delete[] VarDegree;     // 3
    delete[] ChkDegree;     // 4
    dims_2[0] = mxGetM(prhs[5]);
    for (int i=0;i<dims_2[0];i++) delete[] Col_Link_Row[i]; delete[] Col_Link_Row;  // 5
    dims_2[0] = mxGetM(prhs[6]);
    for (int i=0;i<dims_2[0];i++) delete[] Row_Link_Col[i]; delete[] Row_Link_Col;  // 6
    dims_2[0] = mxGetM(prhs[8]);
    for (int i=0;i<dims_2[0];i++) delete[] TABLE_MULTIPLY[i]; delete[] TABLE_MULTIPLY;  // 8
    dims_2[0] = mxGetM(prhs[9]);
    for (int i=0;i<dims_2[0];i++) delete[] TABLE_ADD[i]; delete[] TABLE_ADD;  // 9
    delete[] TABLE_INVERSE;  // 10
    dims_2[0] = mxGetM(prhs[11]);
    for (int i=0;i<dims_2[0];i++) delete[] H[i]; delete[] H;  // 11
    
    dims_3 = (int *)mxGetDimensions(prhs[12]);
    for (int i=0;i<dims_3[0];i++){
        for (int j=0;j<dims_3[1];j++)
            delete[] Q[i][j];
        delete[] Q[i];
    }
    delete[] Q;  //12
    
    dims_3 = (int *)mxGetDimensions(prhs[13]);
    for (int i=0;i<dims_3[0];i++){
        for (int j=0;j<dims_3[1];j++)
            delete[] R[i][j];
        delete[] R[i];
    }
    delete[] R; // 13
    
    dims_2[0] = mxGetM(prhs[14]);
    for (int i=0;i<dims_2[0];i++) delete[] L_Post[i]; delete[] L_Post;
    
    delete[] LDR_Vector;
    delete[] L_SIGMA;
    delete[] L_RHO;
    
    // 输出参数存入plhs
    int dims_2_l[] = {1, CodeLen};
    plhs[0] = mxCreateNumericArray(2, dims_2_l, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(plhs[0]), DecodingSeq, CodeLen*sizeof(unsigned int));
    
    for (int i=0;i<CodeLen;i++)
        Decode_output[i] = (double)DecodingSeq[i];
    
    plhs[1] = mxCreateDoubleMatrix(1, CodeLen, mxREAL);
    ptr_double = (double *)mxGetData(plhs[1]);
    for (int i = 0; i < CodeLen; i++) memcpy(ptr_double+i, Decode_output+i, sizeof(double));
    
    // 输出变量释放空间
    delete[] DecodingSeq;
    delete[] Decode_output;
}

void LLR_BoxPlus(double * L, double * L1, double * L2, unsigned int A1, unsigned int A2,
        unsigned int qAry,
        unsigned int **TABLE_MULTIPLY, unsigned int **TABLE_ADD, unsigned int *TABLE_INVERSE)
{
	double L_temp[256];
	unsigned int alpha_i;
    if (A1 == 0)
    {
        for (unsigned int l = 0; l < qAry - 1; l++)
        {
			alpha_i = l + 1;
            unsigned int first = TABLE_MULTIPLY[TABLE_INVERSE[A2]][alpha_i] - 1;
			L_temp[l]= L2[first];
            
        }
		for (unsigned int l = 0; l < qAry - 1; l++)
        {
			L[l] = L_temp[l];
        }
    }
    else if (A2 == 0)
    {
        for (unsigned int l = 0; l < qAry - 1; l++)
        {
			alpha_i = l + 1;
            unsigned int first = TABLE_MULTIPLY[TABLE_INVERSE[A1]][alpha_i] - 1;
            L_temp[l]= L1[first];
        }
		for (unsigned int l = 0; l < qAry - 1; l++)
        {
			L[l] = L_temp[l];
        }
    }
    else
    {
        double a, b, c, d;
        d = 0;
        for(unsigned int xx = 1; xx < qAry; xx++)
        {
            unsigned int v1 = xx - 1;
            unsigned int v2 = TABLE_MULTIPLY[TABLE_INVERSE[A2]][TABLE_MULTIPLY[xx][A1]] - 1;
            if(d > L1[v1] + L2[v2])
            {
                d = d + log(1 + exp(-1 * (d - ( L1[v1] + L2[v2]) )));
            }
            else
            {
                d = L1[v1] + L2[v2] + log(1 + exp(-1 * (L1[v1] + L2[v2] - d)));
            }
        }
        for (unsigned int l = 0; l < qAry - 1; l++)
        {
			alpha_i = l + 1;
            unsigned int v1 = TABLE_MULTIPLY[TABLE_INVERSE[A1]][alpha_i] - 1;
            unsigned int v2 = TABLE_MULTIPLY[TABLE_INVERSE[A2]][alpha_i] - 1;
            c = 0;
            if(L1[v1] > L2[v2])
            {
                c = L1[v1] + log(1 + exp(-1 * ( L1[v1] - L2[v2] )));
            }
            else
            {
                c = L2[v2] + log(1 + exp(-1 * ( L2[v2] - L1[v1])));
            }
			unsigned int x;
            for (unsigned int xx = 1; xx < qAry; xx++)
            {
				x = xx - 1;
                if (xx != TABLE_MULTIPLY[alpha_i][TABLE_INVERSE[A1]])
                {
                    v1 = xx - 1;
                    v2 = TABLE_MULTIPLY[TABLE_INVERSE[A2]][TABLE_ADD[alpha_i][TABLE_MULTIPLY[xx][A1]]] - 1;
                    if(c > L1[v1] + L2[v2])
                    {
                        c = c + log(1 + exp(-1 * (c - (L1[v1] + L2[v2]))));
                    }
                    else
                    {
                        c =  L1[v1] + L2[v2] + log(1 + exp(-1 * ( L1[v1] + L2[v2] - c )));
                    }
                }
            }
            L_temp[l]= c - d;
        }
		for (unsigned int l = 0; l < qAry - 1; l++)
        {
			L[l] = L_temp[l];
        }
    }
}


bool Decode(double ** LDR, unsigned int *DecodingSeq,    // 0~1
        unsigned long CodeLen, unsigned long ChkLen,      // 2~3
        unsigned int *VarDegree, unsigned int *ChkDegree,   // 4~5
        unsigned long **Col_Link_Row, unsigned long **Row_Link_Col,   // 6~7
        unsigned int qAry, unsigned int **TABLE_MULTIPLY, unsigned int **TABLE_ADD, unsigned int *TABLE_INVERSE,    // 8~11
        unsigned int **H,    // 12
        double ***Q, double  ***R,  // 13~14
        double  **L_Post, double *LDR_Vector, double *L_SIGMA, double *L_RHO,   // 15~18
        unsigned int maxIteration)    // 19
{
    unsigned int j, i, l;
    unsigned int di, dj;
    unsigned int iterations = 0;
    // initial
    // 初始化Q，Q为变量节点的输出消息
    for (i = 0; i < CodeLen; i++)
    {
        for (di = 0; di < VarDegree[i]; di++)
        {
           // mexPrintf("j = Col_Link_Row[%d][%d]\n",i,di);
            j = Col_Link_Row[i][di];
            for (l = 0; l < qAry - 1; l++)
            {
               // mexPrintf("Q[%d][%d][%d] = LDR[%d][%d]\n",j,i,l,i,l);
                // Q[j][i][l] = LDR[i][l];反过来？？
                Q[j][i][l] = LDR[l][i];
            }
        }
    }
    // 初始化R，R为校验节点的输出消息
    for (j = 0; j < ChkLen; j++)
    {
        for (dj = 0; dj < ChkDegree[j]; dj++)
        {
            i = Row_Link_Col[j][dj];
            for (l = 0; l < qAry - 1; l++)
            {
                R[j][i][l] = 0;
            }
        }
    }
    // iterative deocding process
    while(true)
    {
        // tentative coding
        double maxLDR;
        unsigned int alpha_i;
        for (i = 0; i < CodeLen; i ++)
        {
            for (l = 0; l < qAry - 1; l ++)
            {
                //L_Post[i][l] = LDR[i][l];
                L_Post[i][l] = LDR[l][i];
            }
            for (di = 0; di < VarDegree[i]; di++)
            {
                j = Col_Link_Row[i][di];
                for (l = 0; l < qAry - 1; l ++)
                {
                    L_Post[i][l] = L_Post[i][l] + R[j][i][l];
                }
            }
            maxLDR = 0;
            alpha_i = 0;
            for (l = 0; l < qAry - 1; l ++)
            {
                if (maxLDR < L_Post[i][l])
                {
                    maxLDR = L_Post[i][l];
                    alpha_i = l + 1;
                }
            }
            if (maxLDR > 0)
                DecodingSeq[i] = alpha_i;
            else
                DecodingSeq[i] = 0;
			/////////////////////////////////////////////////
		//	if (DecodingSeq[i] != 0)
		//		cout << i + 1 << '\t';
        }
		//cout << '\n';
        // check parity
        unsigned int syndrome = 0;
        bool satisfy_check_parity = true;
        for (j = 0; j < ChkLen; j ++)
        {
            syndrome = 0;
            for (dj = 0; dj < ChkDegree[j]; dj ++)
            {
                i = Row_Link_Col[j][dj];
                syndrome = TABLE_ADD[syndrome][TABLE_MULTIPLY[DecodingSeq[i]][H[j][i]]];
            }
            if (syndrome != 0)
            {
                satisfy_check_parity = false;
                break;
            }
        }
        if (satisfy_check_parity == true)
            return true;
        else if (iterations >= maxIteration)
            return false;
        else
            iterations = iterations + 1;
        // Horizontal Step
        for (i = 0; i < CodeLen; i++)
        {
			for (l = 0; l < qAry - 1; l ++)
            {
				LDR_Vector[l] = 0;
            }
            for (di = 0; di < VarDegree[i]; di ++)
            {
                j = Col_Link_Row[i][di];
                for (l = 0; l < qAry - 1; l ++)
                {
                    LDR_Vector[l] = LDR_Vector[l] + R[j][i][l];
                }
            }
            for (di = 0; di < VarDegree[i]; di ++)
            {
                j = Col_Link_Row[i][di];
                for (l = 0; l < qAry - 1; l ++)
                {
                    //Q[j][i][l] = LDR[i][l] + LDR_Vector[l] - R[j][i][l];
                    Q[j][i][l] = LDR[l][i] + LDR_Vector[l] - R[j][i][l];
                }
            }
        }
        // Vertical Step
        for (j = 0; j < ChkLen; j ++)
        {
            for (dj = 0; dj < ChkDegree[j]; dj ++)
            {
                // forward
                for (l = 0; l < qAry - 1; l ++)
                {
                    L_SIGMA[l] = 0;
                }
                for (unsigned int djj = 0; djj < dj; djj++)
                {
                    i = Row_Link_Col[j][djj];
                    if (djj == 0)
                    {
                        LLR_BoxPlus(L_SIGMA, L_SIGMA, Q[j][i], 0, H[j][i],
                                qAry, TABLE_MULTIPLY, TABLE_ADD, TABLE_INVERSE);
                    }
                    else
                    {
                        LLR_BoxPlus(L_SIGMA, L_SIGMA, Q[j][i], 1, H[j][i],
                                qAry, TABLE_MULTIPLY, TABLE_ADD, TABLE_INVERSE);
                    }
                }
                // backward
                for (l = 0; l < qAry - 1; l ++)
                {
                    L_RHO[l] = 0;
                }
                for (int djj = ChkDegree[j] - 1; djj > dj; djj --)//////////
                {
                    i  = Row_Link_Col[j][djj];
                    if (djj == ChkDegree[j] - 1)
                    {
                        LLR_BoxPlus(L_RHO, L_RHO, Q[j][i], 0, H[j][i],
                                qAry, TABLE_MULTIPLY, TABLE_ADD, TABLE_INVERSE);
                    }
                    else
                    {
                        LLR_BoxPlus(L_RHO, L_RHO, Q[j][i], 1, H[j][i],
                                qAry, TABLE_MULTIPLY, TABLE_ADD, TABLE_INVERSE);
                    }
                }
                // given message
                if (dj == 0)
                {
                    i  = Row_Link_Col[j][dj];
                    LLR_BoxPlus(R[j][i], L_SIGMA, L_RHO, 0, TABLE_INVERSE[H[j][i]],
                            qAry, TABLE_MULTIPLY, TABLE_ADD, TABLE_INVERSE);
                }
                else if (dj == ChkDegree[j] - 1)
                {
                    i  = Row_Link_Col[j][dj];
                    LLR_BoxPlus(R[j][i], L_SIGMA, L_RHO, TABLE_INVERSE[H[j][i]], 0,
                            qAry, TABLE_MULTIPLY, TABLE_ADD, TABLE_INVERSE);
                }
                else
                {
                    i  = Row_Link_Col[j][dj];
                    LLR_BoxPlus(R[j][i], L_SIGMA, L_RHO, TABLE_INVERSE[H[j][i]], TABLE_INVERSE[H[j][i]],
                            qAry, TABLE_MULTIPLY, TABLE_ADD, TABLE_INVERSE);
                }
            }
        }
    }
    return true;
}