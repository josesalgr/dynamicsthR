/*
 * Instancias.cpp
 *
 *  Created on: 08-01-2017
 *      Author: jose
 */

#include "Dynamic_actions.h"
#include <sstream>
#include <string>
#include <limits>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <climits>
#include <limits>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List createModel_DynamicActions(int n_units,
                                int n_species,
                                int n_threats,
                                Rcpp::DataFrame Jk, 
                                Rcpp::DataFrame Ij, 
                                Rcpp::DataFrame Ik,
                                Rcpp::DataFrame ExpansionType,
                                Rcpp::DataFrame Dlong,
                                Rcpp::DataFrame Adyacency,
                                Rcpp::DataFrame Dradial,
                                Rcpp::DataFrame Ck,
                                int levels, 
                                int periods, 
                                int budget_per_period){
  
  std::cout<<"Reading data"<<std::endl;
  
  std::map<int,std::map<int,bool>> matrix_JK;
  
  // Get the number of rows and columns in the data frame
  int numRows = Jk.nrows();
  int numCols = Jk.size();

  // Iterate through rows and columns of the data frame
  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numRows; ++i) {
      // Get the value from the data frame
      Rcpp::IntegerVector df = Jk[j];
      matrix_JK[i][j] = df[i];
    }
  }
  
  std::map<int,std::map<int,bool>> matrix_IJ;
  
  // Get the number of rows and columns in the data frame
  numRows = Ij.nrows();
  numCols = Ij.size();
  
  // Iterate through rows and columns of the data frame
  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numRows; ++i) {
      // Get the value from the data frame
      Rcpp::IntegerVector df = Ij[j];
      matrix_IJ[i][j] = df[i];
    }
  }
  
  std::map<int,std::map<int,int>> matrix_IK;
  
  // Get the number of rows and columns in the data frame
  numRows = Ik.nrows();
  numCols = Ik.size();
  
  // Iterate through rows and columns of the data frame
  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numRows; ++i) {
      // Get the value from the data frame
      Rcpp::IntegerVector df = Ik[j];
      matrix_IK[i][j] = df[i];
    }
  }
  
  std::map<int,std::map<int,int>> matrix_spread;
  
  // Get the number of rows and columns in the data frame
  numRows = ExpansionType.nrows();
  numCols = ExpansionType.size();
  
  // Iterate through rows and columns of the data frame
  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numRows; ++i) {
      // Get the value from the data frame
      Rcpp::NumericVector df = ExpansionType[j];
      matrix_spread[i][j] = df[i];
    }
  }
  
  std::map<int,std::map<int,double>> matrix_Dlong;
  
  // Get the number of rows and columns in the data frame
  numRows = Dlong.nrows();
  numCols = Dlong.size();
  
  // Iterate through rows and columns of the data frame
  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numRows; ++i) {
      // Get the value from the data frame
      Rcpp::NumericVector df = Dlong[j];
      matrix_Dlong[i][j] = df[i];
    }
  }

  std::map<int,std::map<int,double>> matrix_Adj;
  
  // Get the number of rows and columns in the data frame
  numRows = Adyacency.nrows();
  numCols = Adyacency.size();
  
  // Iterate through rows and columns of the data frame
  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numRows; ++i) {
      // Get the value from the data frame
      Rcpp::IntegerVector df = Adyacency[j];
      matrix_Adj[i][j] = df[i];
    }
  }

  std::map<int,std::map<int,double>> matrix_Dradial;
  
  // Get the number of rows and columns in the data frame
  numRows = Dradial.nrows();
  numCols = Dradial.size();
  
  // Iterate through rows and columns of the data frame
  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numRows; ++i) {
      // Get the value from the data frame
      Rcpp::NumericVector df = Dradial[j];
      //std::cout<<"i: "<< i << ", j: "<< j << ", value: " << df[i] <<std::endl;

      matrix_Dradial[i][j] = df[i];
    }
  }

  std::map<int,double> matrix_CK;
  
  // Get the number of rows and columns in the data frame
  numRows = Ck.nrows();
  
  // Iterate through rows and columns of the data frame
  for (int j = 0; j < 1; ++j) {
    for (int i = 0; i < numRows; ++i) {
      // Get the value from the data frame
      Rcpp::NumericVector df = Ck[j];
      matrix_CK[i] = df[i];
    }
  }

  //Parameters
  std::vector<std::string> vtype;
  std::vector<std::string> n_variable;
  std::vector<string> sense;
  std::string name;
  std::vector<double> lb;
  std::vector<double> ub;
  std::vector<double> obj;
  std::vector<double> A_i;
  std::vector<double> A_j;
  std::vector<double> A_x;
  std::vector<double> rhs;
  int row_constraint = 0;
  
  //----------------------------------------------------------------------------
  //Variables definition--------------------------------------------------------
  //----------------------------------------------------------------------------
  
  //X_kitn: Binary variable that indicates whether the actions to abate the threat k with level n 
  //in the unit i in the time t
  for (int i=0;i<n_units;i++)
  {
    for(int k=0;k<n_threats;k++)
    {
      for(int t=0;t<periods;t++)
      {
        for (int n = 0; n < levels; n++)
        {
          vtype.push_back("B");
          name = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n) ;
          n_variable.push_back(name);
          lb.push_back(0);
          ub.push_back(1);
          obj.push_back(0);
        }
      }
    }
  }
  
  //V_kitn: Binary variable that indicates whether exists the threat k with level n 
  //in the unit i in the time t
  for (int i=0;i<n_units;i++)
  {
    for(int k=0;k<n_threats;k++)
    {
      for(int t=0;t<periods;t++)
      {
        for (int n = 0; n < levels; n++)
        {
          vtype.push_back("B");
          name = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n) ;
          n_variable.push_back(name);
          lb.push_back(0);
          ub.push_back(1);
          
          if(n > 0 && t == periods-1){
            obj.push_back(matrix_CK[i]);
          }
          else{
            obj.push_back(0); 
          }
        }
      }
    }
  }
  
  //W_it: Binary variable that indicates if the unit i in the period t is selected to be monitoring
  for (int i=0;i<n_units;i++)
  {
    for (int t = 0; t < periods; t++)
    {
      vtype.push_back("B");
      name = "w_" + std::to_string(i) + "_" + std::to_string(t);
      n_variable.push_back(name);
      lb.push_back(0);
      ub.push_back(1);
      obj.push_back(0);
    }
  }
  
  //Y_ist: Binary variable indicating if the species s exists in the unit i in the period t
  for (int i=0;i<n_units;i++)
  {
    for(int s=0;s<n_species;s++)
    {
      for (int t = 0; t < periods; t++)
      {
        if(matrix_IJ[i][s] > 0){
          vtype.push_back("B");
          name = "y_" + std::to_string(i) + "_" + std::to_string(s) + "_" + std::to_string(t);
          n_variable.push_back(name);
          lb.push_back(0);
          ub.push_back(1);
          
          // if(t == periods-1){
          //   obj.push_back(matrix_CK[i]/n_species);
          // }
          // else{
          obj.push_back(0);
          // }
        }
      }
    }
  }
  
  //Budget per period (variable)
  for (int t=0; t<periods; t++)
  {
    vtype.push_back("C");
    name = "bd_" + std::to_string(t);
    n_variable.push_back(name);
    lb.push_back(0);
    ub.push_back(100000000000);
    obj.push_back(0);
  }
  
  // //Budget default (objetive)
  // vtype.push_back("C");
  // name = "bda";
  // n_variable.push_back(name);
  // lb.push_back(0);
  // ub.push_back(100000000000);
  // obj.push_back(0);
  // 
  // //Min benefit
  // vtype.push_back("C");
  // name = "ben";
  // n_variable.push_back(name);
  // lb.push_back(0);
  // ub.push_back(100000000000);
  // obj.push_back(0);
  
  std::cout<<"Variables created"<<std::endl;
  
  //----------------------------------------------------------------------------
  //Creating initial constraint-------------------------------------------------
  //----------------------------------------------------------------------------
  
  //1) Set the threats in the period 1
  for (int i=0; i<n_units; i++)
  {
    for(int k=0; k<n_threats; k++)
    {
      std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(0) + "_" + std::to_string(0);
      auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
      
      if(matrix_IK[i][k] > 0)
      {
        std::string name_vn = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(0) + "_" + std::to_string(matrix_IK[i][k]);
        auto pos_vn = find(n_variable.begin(), n_variable.end(), name_vn) - n_variable.begin();
        
        A_i.push_back(row_constraint);
        A_j.push_back(pos_vn);
        A_x.push_back(1);
        rhs.push_back(1);
        sense.push_back("==");
      }
      else{
        A_i.push_back(row_constraint);
        A_j.push_back(pos_v);
        A_x.push_back(1);
        rhs.push_back(1);
        sense.push_back("==");
      }
      row_constraint++;
    }
  }
  std::cout<<"Constraint type I created"<<std::endl;
  
  //2) Relation between actions and threats
  for(int i=0; i<n_units; i++)
  {
    for(int k=0; k<n_threats; k++)
    {
      for(int t=0; t<periods; t++)
      {
        for(int n=0; n<levels; n++)
        {
          std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
          std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
          std::string name_v_n0 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(0);
          std::string name_x_n0 = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(0);
          
          auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
          auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
          auto pos_v_n0 = find(n_variable.begin(), n_variable.end(), name_v_n0) - n_variable.begin();
          auto pos_x_n0 = find(n_variable.begin(), n_variable.end(), name_x_n0) - n_variable.begin();
          
          if(n == 0)
          {
            A_i.push_back(row_constraint);
            A_j.push_back(pos_x_n0);
            A_x.push_back(-1);
            
            A_i.push_back(row_constraint);
            A_j.push_back(pos_v_n0);
            A_x.push_back(1);
            
            rhs.push_back(0.0);
            sense.push_back("<=");
            
          }
          else if(n > 0)
          {
            A_i.push_back(row_constraint);
            A_j.push_back(pos_v);
            A_x.push_back(1);
            
            A_i.push_back(row_constraint);
            A_j.push_back(pos_x);
            A_x.push_back(-1);
            
            A_i.push_back(row_constraint);
            A_j.push_back(pos_x_n0);
            A_x.push_back(-1);
            
            rhs.push_back(0.0);
            sense.push_back("<=");
          }
          row_constraint++;
        }
      }
    }
  }
  
  for(int i=0; i<n_units; i++)
  {
    for(int k=0; k<n_threats; k++)
    {
      for(int t=0; t<periods; t++)
      {
        for(int n=0; n<levels; n++)
        {
          std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
          std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
          
          auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
          auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_v);
          A_x.push_back(1);
          
          A_i.push_back(row_constraint + 1);
          A_j.push_back(pos_x);
          A_x.push_back(1);
        }
        rhs.push_back(1.0);
        sense.push_back("==");
        
        rhs.push_back(1.0);
        sense.push_back("==");
        
        row_constraint = row_constraint + 2;
      }
    }
  }
  std::cout<<"Constraint type II created"<<std::endl;
  
  //3) It isn't allow a action in the period 1
  for(int i=0; i<n_units; i++)
  {
    for(int k=0; k<n_threats; k++)
    {
      std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(0) + "_" + std::to_string(0);
      auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
      
      A_i.push_back(row_constraint);
      A_j.push_back(pos_x);
      A_x.push_back(1);
      rhs.push_back(1);
      sense.push_back("==");
      row_constraint++;
    }
  }
  std::cout<<"Constraint type III created"<<std::endl;
  
  //4) Set species distribution at period 1
  for (int i=0; i<n_units; i++)
  {
    for(int s=0; s<n_species; s++)
    {
      if(matrix_IJ[i][s] > 0){
        name = "y_" + std::to_string(i) + "_" + std::to_string(s) + "_" + std::to_string(0);
        auto pos = find(n_variable.begin(), n_variable.end(), name) - n_variable.begin();
        
        A_i.push_back(row_constraint);
        A_j.push_back(pos);
        A_x.push_back(1);
        rhs.push_back(matrix_IJ[i][s]);
        sense.push_back("==");
        row_constraint++;
      }
    }
  }
  std::cout<<"Constraint type IV created"<<std::endl;
  
  //5) Levels
  for (int i=0; i<n_units; i++)
  {
    for(int k=0; k<n_threats; k++)
    {
      for(int t=1; t<periods; t++)
      {
        std::string name_v_nN = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(levels - 1);
        std::string name_v_nN_t1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t - 1) + "_" + std::to_string(levels - 1);
        std::string name_x_nN_t1 = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t - 1) + "_" + std::to_string(levels - 1);
        std::string name_v_n0 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(0);
        std::string name_v_n1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(1);
        
        auto pos_v_nN = find(n_variable.begin(), n_variable.end(), name_v_nN) - n_variable.begin();
        auto pos_v_nN_t1 = find(n_variable.begin(), n_variable.end(), name_v_nN_t1) - n_variable.begin();
        auto pos_x_nN_t1 = find(n_variable.begin(), n_variable.end(), name_x_nN_t1) - n_variable.begin();
        auto pos_v_n0 = find(n_variable.begin(), n_variable.end(), name_v_n0) - n_variable.begin();
        auto pos_v_n1 = find(n_variable.begin(), n_variable.end(), name_v_n1) - n_variable.begin();
        
        for(int n=0; n<levels-1; n++)
        {
          std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t - 1) + "_" + std::to_string(n);
          std::string name_v_np = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n + 1);
          std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t - 1) + "_" + std::to_string(n);
          
          auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
          auto pos_v_np = find(n_variable.begin(), n_variable.end(), name_v_np) - n_variable.begin();
          auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_v);
          A_x.push_back(1);
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_x);
          A_x.push_back(1);
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_v_n0);
          A_x.push_back(-1);
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_v_n1);
          A_x.push_back(-1);
          
          rhs.push_back(1);
          sense.push_back("<=");
          row_constraint++;
          
          if(n>0)
          {
            A_i.push_back(row_constraint);
            A_j.push_back(pos_v);
            A_x.push_back(1);
            
            A_i.push_back(row_constraint);
            A_j.push_back(pos_x);
            A_x.push_back(-1);
            
            A_i.push_back(row_constraint);
            A_j.push_back(pos_v_np);
            A_x.push_back(-1);
            
            rhs.push_back(0);
            sense.push_back("<=");
            row_constraint++;
          }
        }
        
        A_i.push_back(row_constraint);
        A_j.push_back(pos_v_nN_t1);
        A_x.push_back(1);
        
        A_i.push_back(row_constraint);
        A_j.push_back(pos_x_nN_t1);
        A_x.push_back(-1);
        
        A_i.push_back(row_constraint);
        A_j.push_back(pos_v_nN);
        A_x.push_back(-1);
        
        rhs.push_back(0);
        sense.push_back("<=");
        row_constraint++;
        
        for(int n=levels-1; n<=levels-1; n++)
        {
          A_i.push_back(row_constraint);
          A_j.push_back(pos_v_nN_t1);
          A_x.push_back(1);
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_x_nN_t1);
          A_x.push_back(1);
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_v_n1);
          A_x.push_back(-1);
          
          rhs.push_back(1);
          sense.push_back("<=");
          row_constraint++;
        }
      }
    }
  }
  std::cout<<"Constraint type V created"<<std::endl;
  
  //6) Benefits
  for(int k=0; k<n_threats; k++)
  {
    for(int i=0; i<n_units; i++)
    {
      for(int j=0; j<3; j++)
      {
        if(matrix_spread[k][j]>0)
        {
          //Radial propagation
          if(j==0)
          {
            for(int t=0; t<periods-1; t++)
            {
              int Pcounter_units=0;
              std::string name_v_N1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(levels - 1);
              std::string name_v_t1_n1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(1);
              
              auto pos_v_N1 = find(n_variable.begin(), n_variable.end(), name_v_N1) - n_variable.begin();
              auto pos_v_t1_n1 = find(n_variable.begin(), n_variable.end(), name_v_t1_n1) - n_variable.begin();
              
              for(int i_2=0; i_2<n_units; i_2++)
              {
                if((matrix_spread[k][j]>=matrix_Dradial[i][i_2] && matrix_Dradial[i][i_2]>0) || (matrix_Adj[i][i_2]==1))
                {
                  Pcounter_units+=1;
                  
                  std::string name_v_i2_N1 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(levels - 1);
                  std::string name_v_i2_t1_N1 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(levels - 1);
                  
                  auto pos_v_i2_N1 = find(n_variable.begin(), n_variable.end(), name_v_i2_N1) - n_variable.begin();
                  auto pos_v_i2_t1_N1 = find(n_variable.begin(), n_variable.end(), name_v_i2_t1_N1) - n_variable.begin();
                  
                  for(int n=1; n<levels-1; n++)
                  {
                    std::string name_v_i2 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                    std::string name_v_i2_t1 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(n);
                    std::string name_x_i2 = "x_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                    
                    auto pos_v_i2 = find(n_variable.begin(), n_variable.end(), name_v_i2) - n_variable.begin();
                    auto pos_v_i2_t1 = find(n_variable.begin(), n_variable.end(), name_v_i2_t1) - n_variable.begin();
                    auto pos_x_i2 = find(n_variable.begin(), n_variable.end(), name_x_i2) - n_variable.begin();
                    
                    //levels2
                    A_i.push_back(row_constraint);
                    A_j.push_back(pos_v_i2_t1);
                    A_x.push_back(1);
                    
                    //levels3
                    A_i.push_back(row_constraint + 1);
                    A_j.push_back(pos_v_i2);
                    A_x.push_back(-1);
                    
                    A_i.push_back(row_constraint + 1);
                    A_j.push_back(pos_x_i2);
                    A_x.push_back(1);
                    
                  }
                  
                  //levels2
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v_i2_t1_N1);
                  A_x.push_back(1);
                  
                  //levels3
                  A_i.push_back(row_constraint + 1);
                  A_j.push_back(pos_v_i2_N1);
                  A_x.push_back(-1);
                }
              }
              if(Pcounter_units>0)
              {
                for(int n=1; n<levels-1; n++)
                {
                  std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                  std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                  
                  auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
                  auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v);
                  A_x.push_back(-Pcounter_units);
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_x);
                  A_x.push_back(Pcounter_units);
                }
                
                A_i.push_back(row_constraint);
                A_j.push_back(pos_v_N1);
                A_x.push_back(-Pcounter_units);
                
                rhs.push_back(0);
                sense.push_back(">=");
                
                //levels3
                A_i.push_back(row_constraint + 1);
                A_j.push_back(pos_v_t1_n1);
                A_x.push_back(1);
                
                A_i.push_back(row_constraint + 1);
                A_j.push_back(pos_v_N1);
                A_x.push_back(-1);
                
                rhs.push_back(0);
                sense.push_back("<=");
                row_constraint = row_constraint + 2;
              }
            }
          }
          //Downstream propagation
          if(j==1)
          {
            for(int t=0; t<periods-1; t++)
            {
              int Pcounter_units=0;
              int Pcounter_units2=0;
              std::string name_v_N1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(levels - 1);
              std::string name_v_t1_n1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(1);
              std::string name_v_t1_n0 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(0);
              std::string name_v_n0 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(0);
              
              auto pos_v_N1 = find(n_variable.begin(), n_variable.end(), name_v_N1) - n_variable.begin();
              auto pos_v_t1_n1 = find(n_variable.begin(), n_variable.end(), name_v_t1_n1) - n_variable.begin();
              auto pos_v_t1_n0 = find(n_variable.begin(), n_variable.end(), name_v_t1_n0) - n_variable.begin();
              auto pos_v_n0 = find(n_variable.begin(), n_variable.end(), name_v_n0) - n_variable.begin();
              
              for(int i_2=0; i_2<n_units; i_2++)
              {
                std::string name_v_i2_N1 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(levels - 1);
                auto pos_v_i2_N1 = find(n_variable.begin(), n_variable.end(), name_v_i2_N1) - n_variable.begin();
                
                if(((matrix_spread[k][j]>=matrix_Dlong[i][i_2] && matrix_Dlong[i][i_2]!=0)) || (matrix_Adj[i][i_2]==1 && matrix_Dlong[i][i_2]!=0))
                {
                  Pcounter_units+=1;
                  
                  for(int n=1; n<levels; n++)
                  {
                    std::string name_v_i2_t1 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(n);
                    auto pos_v_i2_t1 = find(n_variable.begin(), n_variable.end(), name_v_i2_t1) - n_variable.begin();
                    
                    //levels2
                    A_i.push_back(row_constraint);
                    A_j.push_back(pos_v_i2_t1);
                    A_x.push_back(1);
                  }
                }
                if(((matrix_spread[k][j]>=matrix_Dlong[i_2][i] && matrix_Dlong[i_2][i]!=0)) || (matrix_Adj[i_2][i]==1 && matrix_Dlong[i_2][i]!=0))
                {
                  Pcounter_units2+=1;
                  
                  for(int n=1; n<levels-1; n++)
                  {
                    std::string name_v_i2 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                    std::string name_x_i2 = "x_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                    
                    auto pos_v_i2 = find(n_variable.begin(), n_variable.end(), name_v_i2) - n_variable.begin();
                    auto pos_x_i2 = find(n_variable.begin(), n_variable.end(), name_x_i2) - n_variable.begin();
                    
                    //levels3
                    A_i.push_back(row_constraint + 1);
                    A_j.push_back(pos_v_i2);
                    A_x.push_back(-1);
                    
                    A_i.push_back(row_constraint + 1);
                    A_j.push_back(pos_x_i2);
                    A_x.push_back(1);
                  }
                  A_i.push_back(row_constraint + 1);
                  A_j.push_back(pos_v_i2_N1);
                  A_x.push_back(-1);
                }
              }
              if(Pcounter_units>0)
              {
                for(int n=1; n<levels-1; n++)
                {
                  std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                  std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                  
                  auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
                  auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v);
                  A_x.push_back(-Pcounter_units);
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_x);
                  A_x.push_back(Pcounter_units);
                }
                A_i.push_back(row_constraint);
                A_j.push_back(pos_v_N1);
                A_x.push_back(-Pcounter_units);
                
                rhs.push_back(0);
                sense.push_back(">=");
              }
              if(Pcounter_units2>0)
              {
                //levels3
                A_i.push_back(row_constraint + 1);
                A_j.push_back(pos_v_t1_n1);
                A_x.push_back(1);
                
                A_i.push_back(row_constraint + 1);
                A_j.push_back(pos_v_N1);
                A_x.push_back(-1);
                
                rhs.push_back(0);
                sense.push_back("<=");
                row_constraint = row_constraint + 2;
              }
              else if(Pcounter_units2==0 && matrix_IK[i][k]==0)
              {
                if(t==periods-2)
                {
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v_t1_n0);
                  A_x.push_back(1);
                  rhs.push_back(1);
                  sense.push_back("==");
                  row_constraint++;
                }
                
                A_i.push_back(row_constraint);
                A_j.push_back(pos_v_n0);
                A_x.push_back(1);
                rhs.push_back(1);
                sense.push_back("==");
                row_constraint++;
              }
              else if(Pcounter_units2==0 && matrix_IK[i][k]>0)
              {
                for(int n=1; n<levels; n++)
                {
                  std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                  std::string name_v_t1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(n);
                  
                  auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
                  auto pos_v_t1 = find(n_variable.begin(), n_variable.end(), name_v_t1) - n_variable.begin();
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v_t1);
                  A_x.push_back(1);
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v);
                  A_x.push_back(-1);
                }
                rhs.push_back(0);
                sense.push_back("<=");
                row_constraint++;
              }
            }
          }
          //Upstream propagation
          if(j==2)
          {
            for(int t=0; t<periods-1; t++)
            {
              int Pcounter_units=0;
              int Pcounter_units2=0;
              std::string name_v_N1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(levels - 1);
              std::string name_v_t1_n1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(1);
              std::string name_v_t1_n0 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(0);
              std::string name_v_n0 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(0);
              
              auto pos_v_N1 = find(n_variable.begin(), n_variable.end(), name_v_N1) - n_variable.begin();
              auto pos_v_t1_n1 = find(n_variable.begin(), n_variable.end(), name_v_t1_n1) - n_variable.begin();
              auto pos_v_t1_n0 = find(n_variable.begin(), n_variable.end(), name_v_t1_n0) - n_variable.begin();
              auto pos_v_n0 = find(n_variable.begin(), n_variable.end(), name_v_n0) - n_variable.begin();
              
              for(int i_2=0; i_2<n_units; i_2++)
              {
                std::string name_v_i2_N1 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(levels - 1);
                auto pos_v_i2_N1 = find(n_variable.begin(), n_variable.end(), name_v_i2_N1) - n_variable.begin();
                
                if(((matrix_spread[k][j]>=matrix_Dlong[i_2][i] && matrix_Dlong[i_2][i]!=0)) || (matrix_Adj[i_2][i]==1 && matrix_Dlong[i_2][i]!=0))
                {
                  Pcounter_units+=1;
                  
                  for(int n=1; n<levels; n++)
                  {
                    std::string name_v_i2_t1 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(n);
                    auto pos_v_i2_t1 = find(n_variable.begin(), n_variable.end(), name_v_i2_t1) - n_variable.begin();
                    
                    //levels2
                    A_i.push_back(row_constraint);
                    A_j.push_back(pos_v_i2_t1);
                    A_x.push_back(1);
                  }
                }
                if(((matrix_spread[k][j]>=matrix_Dlong[i][i_2] && matrix_Dlong[i][i_2]!=0)) || (matrix_Adj[i][i_2]==1 && matrix_Dlong[i][i_2]!=0))
                {
                  Pcounter_units2+=1;
                  
                  for(int n=1; n<levels-1; n++)
                  {
                    std::string name_v_i2 = "v_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                    std::string name_x_i2 = "x_" + std::to_string(i_2) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                    
                    auto pos_v_i2 = find(n_variable.begin(), n_variable.end(), name_v_i2) - n_variable.begin();
                    auto pos_x_i2 = find(n_variable.begin(), n_variable.end(), name_x_i2) - n_variable.begin();
                    
                    //levels3
                    A_i.push_back(row_constraint + 1);
                    A_j.push_back(pos_v_i2);
                    A_x.push_back(-1);
                    
                    A_i.push_back(row_constraint + 1);
                    A_j.push_back(pos_x_i2);
                    A_x.push_back(1);
                  }
                  A_i.push_back(row_constraint + 1);
                  A_j.push_back(pos_v_i2_N1);
                  A_x.push_back(-1);
                }
              }
              if(Pcounter_units>0)
              {
                for(int n=1; n<levels-1; n++)
                {
                  std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                  std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                  
                  auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
                  auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v);
                  A_x.push_back(-Pcounter_units);
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_x);
                  A_x.push_back(Pcounter_units);
                }
                A_i.push_back(row_constraint);
                A_j.push_back(pos_v_N1);
                A_x.push_back(-Pcounter_units);
                
                rhs.push_back(0);
                sense.push_back(">=");
              }
              if(Pcounter_units2>0)
              {
                //levels3
                A_i.push_back(row_constraint + 1);
                A_j.push_back(pos_v_t1_n1);
                A_x.push_back(1);
                
                A_i.push_back(row_constraint + 1);
                A_j.push_back(pos_v_N1);
                A_x.push_back(-1);
                
                rhs.push_back(0);
                sense.push_back("<=");
                row_constraint = row_constraint + 2;
              }
              else if(Pcounter_units2==0 && matrix_IK[i][k]==0)
              {
                if(t==periods-2)
                {
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v_t1_n0);
                  A_x.push_back(1);
                  rhs.push_back(1);
                  sense.push_back("==");
                  row_constraint++;
                }
                
                A_i.push_back(row_constraint);
                A_j.push_back(pos_v_n0);
                A_x.push_back(1);
                rhs.push_back(1);
                sense.push_back("==");
                row_constraint++;
              }
              else if(Pcounter_units2==0 && matrix_IK[i][k]>0)
              {
                for(int n=1; n<levels; n++)
                {
                  std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
                  std::string name_v_t1 = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t + 1) + "_" + std::to_string(n);
                  
                  auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
                  auto pos_v_t1 = find(n_variable.begin(), n_variable.end(), name_v_t1) - n_variable.begin();
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v_t1);
                  A_x.push_back(1);
                  
                  A_i.push_back(row_constraint);
                  A_j.push_back(pos_v);
                  A_x.push_back(-1);
                }
                rhs.push_back(0);
                sense.push_back("<=");
                row_constraint++;
              }
            }
          }
        }
      }
    }
  }
  
  std::cout<<"Constraint type VI created"<<std::endl;
  
  
  //Elimination of species when coexist with threats
  for(int i=0; i<n_units; i++)
  {
    for(int k=0; k<n_threats; k++)
    {
      for(int t=1; t<periods; t++)
      {
        int Psum_species=0;
        
        for(int s=0; s<n_species; s++)
        {
          std::string name_y = "y_" + std::to_string(i) + "_" + std::to_string(s) + "_" + std::to_string(t);
          auto pos_y = find(n_variable.begin(), n_variable.end(), name_y) - n_variable.begin();
          
          if(matrix_IJ[i][s]>0)
          {
            A_i.push_back(row_constraint);
            A_j.push_back(pos_y);
            A_x.push_back(1);
            Psum_species+=1;
          }
        }
        if(Psum_species>0)
        {
          for(int n=1; n<levels; n++)
          {
            std::string name_v = "v_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t - 1) + "_" + std::to_string(n);
            std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t - 1) + "_" + std::to_string(n);
            
            auto pos_v = find(n_variable.begin(), n_variable.end(), name_v) - n_variable.begin();
            auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
            
            A_i.push_back(row_constraint);
            A_j.push_back(pos_v);
            A_x.push_back(Psum_species);
            
            A_i.push_back(row_constraint);
            A_j.push_back(pos_x);
            A_x.push_back(-Psum_species);
            
            rhs.push_back(Psum_species);
            sense.push_back("<=");
            row_constraint++;
          }
        }
      }
    }
  }
  
  std::cout<<"Constraint type VII created"<<std::endl;
  
  //FALTA AGREGAR RECOLONIZACIÓN AQUI
  
  
  
  //Do monitoring after apply actions in a specific unit
  for(int i=0; i<n_units; i++)
  {
    for(int t=1; t<periods; t++)
    {
      std::string name_w = "w_" + std::to_string(i) + "_" + std::to_string(t - 1);
      
      auto pos_w = find(n_variable.begin(), n_variable.end(), name_w) - n_variable.begin();
      
      for(int k=0; k<n_threats; k++)
      {
        for(int n=1; n<levels; n++)
        {
          std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
          
          auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_x);
          A_x.push_back(1);
        }
      }
      A_i.push_back(row_constraint);
      A_j.push_back(pos_w);
      A_x.push_back(-n_threats);
      
      rhs.push_back(0);
      sense.push_back("<=");
      row_constraint++;
    }
  }
  
  //Do monitoring before apply actions in a specific unit
  for(int i=0; i<n_units; i++)
  {
    for(int t=1; t<periods-1; t++)
    {
      std::string name_w = "w_" + std::to_string(i) + "_" + std::to_string(t + 1);
      
      auto pos_w = find(n_variable.begin(), n_variable.end(), name_w) - n_variable.begin();
      
      for(int k=0; k<n_threats; k++)
      {
        for(int n=1; n<levels; n++)
        {
          std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
          
          auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_x);
          A_x.push_back(1);
        }
      }
      A_i.push_back(row_constraint);
      A_j.push_back(pos_w);
      A_x.push_back(-n_threats);
      
      rhs.push_back(0);
      sense.push_back("<=");
      row_constraint++;
    }
  }
  
  std::cout<<"Constraint type VIII created"<<std::endl;
  
  
  std::string name_bd = "bd_" + std::to_string(0);
  auto pos_bd = find(n_variable.begin(), n_variable.end(), name_bd) - n_variable.begin();
  
  A_i.push_back(row_constraint);
  A_j.push_back(pos_bd);
  A_x.push_back(1);
  rhs.push_back(budget_per_period);
  sense.push_back("==");
  row_constraint++;
  
  for(int t=0; t<periods; t++)
  {
    for(int i=0; i<n_units; i++)
    {
      std::string name_w = "w_" + std::to_string(i) + "_" + std::to_string(t);
      auto pos_w = find(n_variable.begin(), n_variable.end(), name_w) - n_variable.begin();
      
      A_i.push_back(row_constraint);
      A_j.push_back(pos_w);
      A_x.push_back(matrix_CK[i]);
      
      for(int k=0; k<n_threats; k++)
      {
        for(int n=1; n<levels; n++)
        {
          std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
          auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
          
          A_i.push_back(row_constraint);
          A_j.push_back(pos_x);
          A_x.push_back(n*matrix_CK[i]);
        }
      }
    }
    
    std::string name_bd = "bd_" + std::to_string(t);
    auto pos_bd = find(n_variable.begin(), n_variable.end(), name_bd) - n_variable.begin();
    
    A_i.push_back(row_constraint);
    A_j.push_back(pos_bd);
    A_x.push_back(-1);
    rhs.push_back(0);
    sense.push_back("<=");
    row_constraint++;
  }
  
  for(int t=0; t<periods; t++)
  {
    if(t<periods-1)
    {
      for(int i=0; i<n_units; i++)
      {
        std::string name_w = "w_" + std::to_string(i) + "_" + std::to_string(t);
        auto pos_w = find(n_variable.begin(), n_variable.end(), name_w) - n_variable.begin();
        
        A_i.push_back(row_constraint);
        A_j.push_back(pos_w);
        A_x.push_back(matrix_CK[i]);
        
        for(int k=0; k<n_threats; k++)
        {
          for(int n=1; n<levels; n++)
          {
            std::string name_x = "x_" + std::to_string(i) + "_" + std::to_string(k) + "_" + std::to_string(t) + "_" + std::to_string(n);
            auto pos_x = find(n_variable.begin(), n_variable.end(), name_x) - n_variable.begin();
            
            A_i.push_back(row_constraint);
            A_j.push_back(pos_x);
            A_x.push_back(n*matrix_CK[i]);
          }
        }
      }
      
      std::string name_bd = "bd_" + std::to_string(t);
      std::string name_bd_t1 = "bd_" + std::to_string(t + 1);
      auto pos_bd = find(n_variable.begin(), n_variable.end(), name_bd) - n_variable.begin();
      auto pos_bd_t1 = find(n_variable.begin(), n_variable.end(), name_bd_t1) - n_variable.begin();
      
      A_i.push_back(row_constraint);
      A_j.push_back(pos_bd_t1);
      A_x.push_back(1);
      
      A_i.push_back(row_constraint);
      A_j.push_back(pos_bd);
      A_x.push_back(-1);
      
      rhs.push_back(budget_per_period);
      sense.push_back("==");
      row_constraint++;
    }
  }
  
  std::cout<<"Constraint type IX created"<<std::endl;
  
  
  //FALTA AGREGAR TEMPORAL CONSTRAINT AQUI
  
  
  return List::create(Named("obj") = obj,
                      Named("n_variable") = n_variable,
                      Named("rhs") = rhs,
                      Named("sense") = sense,
                      Named("A_i") = A_i,
                      Named("A_j") = A_j,
                      Named("A_x") = A_x,
                      Named("modelsense") = "min",
                      Named("lb") = lb,
                      Named("ub") = ub,
                      Named("vtype") = vtype);
}


// Register functions
RCPP_MODULE(yada) {
  Rcpp::function("readProblem" , &createModel_DynamicActions);
}
