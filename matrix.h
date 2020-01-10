/*
    matrix.h
    copyright Purelight. Chan Zee Lok.
*/

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <math.h>
using namespace std;
class matrix {
    public:
        /*initialize*/
        matrix(int,int);
        matrix();
        matrix(int,int,double **);
        matrix(const matrix &);
        matrix & operator=(const matrix &);
        /*iostream*/
        friend istream &operator>>(istream &,matrix &);
        friend ostream &operator<<(ostream &,matrix &);
        /*operator*/
        friend matrix operator+(const matrix &,const matrix &);
        friend matrix operator-(const matrix &,const matrix &);
        friend matrix operator*(const matrix &,const matrix &);
        friend matrix operator*(const matrix &,const double &);
        friend matrix operator*(const double &,const matrix &);
        friend matrix operator/(const matrix &,const double &);
        /*self_operation*/
        matrix & operator+=(const matrix &);
        matrix & operator-=(const matrix &);
        matrix & operator*=(const matrix &);
        matrix & operator*=(const double);
        matrix & operator/=(const double);
        matrix  operator^(const int);
        /*basic operation*/
        void swap_rows(const int,const int);
        void row_product_a_const(const int,const double);
        void add_line_i_to_line_j(const int,const int,const double=1);
        void swap_list(const int,const int);
        void list_product_a_const(const int,const double);
        void add_list_i_to_list_j(const int,const int,const double=1);
        matrix transpose();
        /*advanced opreation*/
        double cofactor(const int,const int);//余子式
        double algebraic_cofactor(const int,const int); //代数余子式
        matrix echelon();//阶梯型
        matrix inverse(); //矩阵的逆
        /*Nature of the matrix*/
        double det();//行列式
        int rank();//秩
        double * operator[](const int);
    private:
        int n,m; //size of the matrix
        int r; //rank of the matrix
        double ** p;
        void alloc();
        void clear();
        void calculate_rank();
        matrix calculate_multiplys(const matrix base,const int a);
};

/*function for alloc space and free space*/
void matrix::alloc()
{
    p=new double*[n];
    for (int i=0; i<n; i++)
        p[i]=new double[m];
}
void matrix::clear()
{
    for (int i=0; i<n; i++)
        if (p[i]!=NULL)
            delete []p[i];
    if (p!=NULL)
        delete []p;
}

/*function for initialize the matrix*/
matrix::matrix()
{
    p=NULL;
    n=m=r=0;
}
matrix::matrix(int a,int b):
    n(a),m(b),r(0),p(NULL)
{
    alloc();
}
matrix::matrix(int a,int b,double ** c):
    n(a),m(b)
{
    alloc();
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            p[i][j]=c[i][j];
}
matrix::matrix(const matrix & a):
    n(a.n),m(a.m)
{
    alloc();
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            p[i][j]=a.p[i][j];
}
matrix & matrix::operator=(const matrix & a)
{
    clear();
    n=a.n;
    m=a.m;
    r=a.r;
    alloc();
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            p[i][j]=a.p[i][j];
    return *this;
}

/*function for iostream*/
istream &operator>>(istream &input,matrix &a)
{
    for (int i=0;i<a.n;i++)
        for (int j=0;j<a.m;j++)
            input>>a[i][j];
    return input;
}
ostream &operator<<(ostream &output,matrix &a)
{
    for (int i=0;i<a.n;i++){
        for (int j=0;j<a.m;j++)
            output<<a.p[i][j]<<' ';
        output<<endl;
    }
    return output;
}

/*function for self-operation of matrix*/
matrix & matrix::operator+=(const matrix & a)
{
    if (n!=a.n||m!=a.m) {
        cerr<<"Error! The two matrix must have the same size when add together!"<<endl;
        exit(0);
    }
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            p[i][j]+=a.p[i][j];
    return *this;
}
matrix & matrix::operator-=(const matrix & a)
{
    if (n!=a.n||m!=a.m) {
        cerr<<"Error! The two matrix must have the same size when subtrace together!"<<endl;
        exit(0);
    }
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            p[i][j]-=a.p[i][j];
    return *this;
}
matrix & matrix::operator*=(const matrix & a)
{
    if (m!=a.n) {
        cerr<<"Error! The number of list in matrix 1 must be the same as the number of rows in matrix 2 when they product!"<<endl;
        exit(0);
    }
    matrix temp(n,a.m);
    for (int i=0; i<n; i++)
        for (int j=0; j<a.m; j++) {
            temp.p[i][j]=0;
            for (int k=0; k<m; k++)
                temp.p[i][j]+=p[i][k]*a.p[k][j];
        }
    return (*this=temp);
}
matrix & matrix::operator*=(const double a)
{
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            p[i][j]*=a;
    return *this;
}
matrix & matrix::operator/=(const double a)
{
    if (a==0) {
        cerr<<"0 cannot be the divisor!"<<endl;
        exit(0);
    }
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            p[i][j]/=a;
    return *this;
}
matrix matrix::operator^(const int a)
{
    if (a<=0) {
        cerr<<"When you use pow in matrix, please be aware that the exp should bigger than 0!"<<endl;
        exit(0);
    }
    return calculate_multiplys(*this,a);
}

/*function for basic operation of matrix*/
matrix operator+(const matrix & a,const matrix & b)
{
    matrix temp(a);
    return (temp+=b);
}
matrix operator-(const matrix & a,const matrix & b)
{
    matrix temp(a);
    return (temp-=b);
}
matrix operator*(const matrix & a,const matrix & b)
{
    matrix temp(a);
    return (temp*=b);
}
matrix operator*(const matrix & a,const double & b)
{
    matrix temp(a);
    return (temp*=b);
}
matrix operator*(const double & a,const matrix & b)
{
    matrix temp(b);
    return (temp*=a);
}
matrix operator/(const matrix & a,const double & b)
{
    matrix temp(a);
    return (temp/=b);
}

/*basic operation of a matrix*/
void matrix::swap_rows(const int a,const int b) //swap ath & bth rows
{
    if (a>n||b>n||a<=0||b<=0) {
        cerr<<"Error! swap_rows(int,int) exceed! "<<endl;
        exit(0);
    }
    for (int k=0; k<m; k++) {
        double temp=p[a-1][k];
        p[a-1][k]=p[b-1][k];
        p[b-1][k]=temp;
    }
}
void matrix::row_product_a_const(const int a,const double b) //ath row * b
{
    if (a>n||a<=0) {
        cerr<<"Error! row_product_a_const(int,double) exceed!"<<endl;
        exit(0);
    }
    for (int k=0; k<m; k++) {
        p[a-1][k]*=b;
        if (p[a-1][k]==-0)
            p[a-1][k]=0;
    }
}
void matrix::add_line_i_to_line_j(const int a,const int b,const double c)
{
    if (a>n||b>n||a<=0||b<=0) {
        cerr<<"Error! add_line_i_to_line_j(int,int,double) exceed! "<<endl;
        exit(0);
    }
    for (int k=0; k<m; k++) {
        p[b-1][k]+=c*p[a-1][k];
        if (p[b-1][k]==-0)
            p[b-1][k]=0;
    }
}
void matrix::swap_list(const int a,const int b) //swap ath & bth list
{
    if (a>m||b>m||a<=0||b<=0) {
        cerr<<"Error! swap_list(int,int) exceed! "<<endl;
        exit(0);
    }
    for (int k=0; k<n; k++) {
        double temp=p[k][a-1];
        p[k][a-1]=p[k][b-1];
        p[k][b-1]=temp;
    }
}
void matrix::list_product_a_const(const int a,const double b) //ath row * b
{
    if (a>m||a<=0) {
        cerr<<"Error! list_product_a_const(int,double) exceed!"<<endl;
        exit(0);
    }
    for (int k=0; k<n; k++) {
        p[k][a-1]*=b;
        if (p[k][a-1]==-0)
            p[k][a-1]=0;
    }
}
void matrix::add_list_i_to_list_j(const int a,const int b,const double c) //j=j+c*i
{
    if (a>m||b>m||a<=0||b<=0) {
        cerr<<"Error! add_list_i_to_list_j(int,int,double) exceed! "<<endl;
        exit(0);
    }
    for (int k=0; k<n; k++) {
        p[k][b-1]+=c*p[k][a-1];
        if (p[k][b-1]==-0)
            p[k][b-1]=0;
    }
}
matrix matrix::transpose()
{
    matrix temp(m,n);
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            temp.p[j][i]=p[i][j];
    temp.r=r;
    return temp;
}

/*advanced operations*/
double matrix::cofactor(const int a,const int b)
{
    if (a>m||b>m||a<=0||b<=0) {
        cerr<<"Error! cofactor(int,int,double) exceed! "<<endl;
        exit(0);
    }
    matrix temp(n-1,m-1);
    for (int i=0,k=0; i<n; i++)
        if (i!=a-1) {
            for (int j=0,l=0; j<m; j++)
                if (j!=b-1)
                    temp.p[k][l++]=p[i][j];
            k++;
        }
    return temp.det();
}
double matrix::algebraic_cofactor(const int a,const int b)
{
    if (a>m||b>m||a<=0||b<=0) {
        cerr<<"Error! algebraic_cofactor(int,int,double) exceed! "<<endl;
        exit(0);
    }
    if ((a+b)%2)
        return -1*cofactor(a,b);
    else
        return cofactor(a,b);
}
matrix matrix::echelon()
{
    matrix temp=*this;
    for (int j=0; j<m-1&&j<n; j++) { //process jth list
        if (temp.p[j][j]==0) {
            for (int i=j+1; i<n; i++)
                if (temp.p[i][j]!=0) {
                    temp.swap_rows(i+1,j+1);
                    break;
                }
        }
        for (int i=j+1; i<n; i++)
            if (temp.p[i][j]!=0) {
                double c=-1*temp.p[i][j]/temp.p[j][j];
                temp.add_line_i_to_line_j(j+1,i+1,c);
            }
    }
    for (int i=0; i<n; i++) {
        int j=i;
        while (temp.p[i][j]==0)
            j++;
        temp.row_product_a_const(i+1,1.0/temp.p[i][j]);
    }
    return temp;
}
matrix matrix::inverse()
{
    if (n!=m||rank()!=n) {
        cerr<<"No inverse matrix!"<<endl;
        exit(0);
    }
    matrix temp(n,m);
    for (int i=0;i<n;i++)
        for (int j=0;j<m;j++)
            temp[i][j]=algebraic_cofactor(j+1,i+1);
    return temp/det();
}

/*Nature of the matrix*/
double matrix::det()
{
    if (n!=m) {
        cerr<<"The matrix has no det!"<<endl;
        exit(0);
    }
    if (n==1)
        return p[0][0];
    double ans=0;
    for (int i=0; i<m; i++)
        ans+=p[0][i]*algebraic_cofactor(1,i+1);
    return ans;
}
int matrix::rank()
{
    calculate_rank();
    return r;
}
double * matrix::operator[](const int a)
{
    return p[a];
}
void matrix::calculate_rank()
{
    matrix temp=echelon();
    if (n<=0||m<=0)
        r=0;
    r=0;
    int pre=-1,now=0;
    for (int i=0;i<n;i++){
        while (now<=m&&temp[i][now]==0)
            now++;
        if (now!=pre){
            pre=now;
            r++;
        }
        now=0;
    }
}
matrix matrix::calculate_multiplys(const matrix base,const int a)
{
    if (a==1)
        return base;
	if (a == 2)
		return base * base;
    if (a%2==0)
        return calculate_multiplys(base,a/2)^2;
    return base*calculate_multiplys(base,a/2)^2;
}

#endif
