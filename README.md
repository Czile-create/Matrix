# Matrix
A C++ class for matrix, containing function such as inverse,det,rank,echelon and several operations for matrix.
## How to use?
### 1. Define a matrix: 
` matrix p(n,m); `  
n,m represent the number of line(somebody may call it as 'row') and list.
### 2. To read/print a matrix: 
`#include <iostream>
 std::cin>>p;
 std::cout<<p;`
### 3. To do some operations: 
  `p1+p2,p1-p2,p1*p2,p1*const,p1/const,p1^const`  
  All these operation can be used.
### 4. To do some advanced operations:
  Here are a list of operations you can use:  
  
  `
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
    ` 
        
## LICENSE
  Please see the file "License"
## Contributor
  Czile copyright.2020
