/*---------------------------------------------------------------------------------------------------------/
g++ -g -o `root-config --cflags` `lhapdf-config --cflags` `gsl-config --cflags` res.cpp `root-config --glibs` `lhapdf-config --libs` `gsl-config --libs` -I ./ -lstdc++ -lm -lMathMore -lLHAPDF -o res && ./res

 //mode: -1 default, 0 adaptative, 1 vegas(MC), 2 miser(MC), 3 plain(MC)
/---------------------------------------------------------------------------------------------------------*/


#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TTree.h"
#include <TFile.h>

#include "Math/AdaptiveIntegratorMultiDim.h"
#include <Math/WrappedParamFunction.h>

#include "Math/IntegratorMultiDim.h"
#include "Math/Functor.h"
#include "Math/IFunctionfwd.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/Integrator.h"
#include <Math/IFunction.h>
#include <Math/GSLMCIntegrator.h>

#include <cmath>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/GridPDF.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include <TMath.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TApplication.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <TRandom3.h>
#include <vector>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <iostream>
#include <random>
#include <cmath>
#include "TCanvas.h"
#include "TGraph.h"
#include <iostream>
#include <TMath.h>
#include <TRandom3.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGaxis.h>



using namespace LHAPDF;
using namespace std;

     //   PDF  //
PDF* nnff = mkPDF("NNPDF21_lo_as_0119_100");


/*----------------------------------------------------------------------------------------------------*/

// Função para obter a integral usando o pacote de integração do ROOT

double GetIntegralByRoot(double (*funcao)(const double *,const double *), int fdim, double *parr, int pdim, double *xa, double *xb, int intmode){

	// kADAPTIVE, kVEGAS, kPLAIN, kMISER
	ROOT::Math::WrappedParamFunction<> integrando(*funcao, fdim, parr, parr+pdim);
	ROOT::Math::AdaptiveIntegratorMultiDim integral(intmode);
	integral.SetFunction(integrando);
	//integral.SetAbsTolerance(1e-3);
	//integral.SetRelTolerance(1e-9);
	//integral.SetSize(1e9);
	integral.SetMaxPts(1e9);
	integral.Integral(xa, xb);

	return integral.Result();
}



double GetIntegralByRoot(double (*funcao)(const double *,const double *), int fdim, double *parr, int pdim, double *xa, double *xb){
	return GetIntegralByRoot(*funcao,fdim,parr,pdim,xa,xb);
}


/*----------------------------------------------------------------------------------------------------*/
                                          /*EXEMPLOS*/
 
/* -----------------------------------------------------------------------------------------------------


// // Função a ser integrada - GSL : f(x, y, z) = 2x * y^2 * e^x
double my_function_GSL(double x[], size_t dim, void* params) {
    return 2.0 * x[0] * std::pow(x[1], 2) * std::exp(x[2]);;
}




// Função a ser integrada para ROOT  X * Y
double integg(const double *var , const double *p){
	double x = var[0];
	double y = var[1];
	double p1 = p[0];
	
	return x*y;
	

}


// Função a ser integrada ROOT 
double my_function_ROOT(const double* x, const double* p) {
    double x_val = x[0]; // Valor da primeira variável x
    double y_val = x[1]; // Valor da segunda variável y
    double z_val = x[2]; // Valor da terceira variável z

    return 2.0 * x_val * y_val * y_val * std::exp(z_val);
}



*//*----------------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------*/




// Função da seção de choque diferencial -ROOT
double differential_cross_section(const double *vars, const double *params) {

    double y3 = vars[0];
    double y4 = vars[1];
    

    // double x1fi = nnff->xfxQ(21,x1,mu2);
   // double x2fj = nnff->xfxQ(21,x2,mu2);
    
    // Parâmetros físicos
 
    
    double m = params[0];
    double pT2 = params[1];
    
    double mT = std::sqrt( pT2 * pT2  + m * m);
    double Mu2 = std::pow(mT,4);
 
///////////////////////////////////////////////////////////////////   
    
    // Cálculo de g usando as relações fornecidas
   
    int nf = 5; // Número de sabores 
    double beta = 11.0 - (2.0 / 3.0) * nf;
    
    //lambda
    double lambda = 0.2;
    
    // Cálculo de alpha_s(Q^2)
    double alpha_s = 4 * M_PI / (beta * log( std::pow(mT,2) / std::pow(lambda,2)));
    
   
    double g2 = 4 * M_PI * alpha_s;
    double g = sqrt(g2);
    
//////////////////////////////////////////////////////////////////
    
    // X1 e x2
    
   // double s = sqrt ( 2 * std::pow(mT,2) * (1 + std::cosh(y3 - y4)) );
    double s = 10000 ;
    double x1 = (mT / sqrt(s))  * (exp(y3) + exp(y4));
    double x2 = (mT / sqrt(s))  * (exp(-y3) + exp(-y4));
    
    
    // Expressão da seção de choque diferencial
    double term1 = 1.0 / (64 * M_PI * std::pow(mT, 4) * std::pow(1 + std::cosh(y3 - y4), 2));
    double term2 = (std::pow(g, 4) / 24) * ((8 * std::cosh(y3 - y4) - 1) / (1 + std::cosh(y3 - y4)));
    double term3 = (std::cosh(y3 - y4) + 2 * std::pow(m, 2) / std::pow(mT, 2) - 2 * std::pow(m, 4) / std::pow(mT, 4));
    double term4 = (4 * std::pow(g, 4) / 9) * (1 / (1 + std::cosh(y3 - y4))) * (std::cosh(y3 - y4) + std::pow(m, 2) / std::pow(mT, 2));
    
    


    
    
    double n = 0;
    double term5= 0;
    double Sterm5=0;

    for (n = 0; n < 7; n += 1) {
        double fi = nnff->xfxQ( n, x1, Mu2);
        double fj = nnff->xfxQ(-n, x2, Mu2);
        
       

        Sterm5 = fj * x2 * fi * x1;
        std::cout << "X1 = " <<  x1 << std::endl;
        term5 = (Sterm5 + term5); 
       // std::cout << "termo 5= " << term5 << std::endl;
    // Faça algo com term5, se necessário.
}    
    
    
    
    return term1 * term2 * term3 * term4 *term5 ;
}



/*

g^{2}=4\pi \alpha_{s}  ---->  \mu= M_{t}^{2} \ \ ||| \ \ \beta = 11- \frac{2}{3} n_{f} ----->\alpha_{s} (Q^{2}) = \frac{4 \pi}{\beta ln(\frac{Q^{2}}{\Lambda ^{2}})} ----> \Lambda ^{2} = \mu^{2} e^{- \frac{4\pi}{\beta \alpha({\mu^{2})}}}

*/



/*
// Função que descreve a seção de choque diferencial
double differential_cross_section_GSL(double* x, size_t dim, void* params) {
    // Defina suas constantes e variáveis aqui
    double g = 1.0; // Exemplo: constante de acoplamento
    double m = 1.0; // Exemplo: massa da partícula
    double y3 = x[0];
    double y4 = x[1];
    double pt = 100; // Convertendo PT^2 para PT

    // Cálculos intermediários da expressão
    double mT = sqrt(2 + m * m);
    double numerator = (cosh(y3 - y4) + 2.0 * m * m / (mT * mT) - 2.0 * m * m * m * m / (mT * mT * mT * mT));
    double denominator = 64.0 * M_PI * mT * mT * mT * mT * pow((1.0 + cosh(y3 - y4)), 2) * (cosh(y3 - y4) + 2.0 * m * m / (mT * mT));

    // Função de espalhamento
    double scattering_function = g * g * g * g / 24.0 * (8.0 * cosh(y3 - y4) - 1.0) / (1.0 + cosh(y3 - y4)) * numerator / denominator;

    return scattering_function;
}

*/


int main(){  
    
    
    /* --------------------------- GSL-------------------------------*/
   
   /*
    const gsl_rng_type* T;
    gsl_rng* r;

    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_function G = { &my_function_GSL, 3, 0 };
    
    double result, error;

    // Limites de integração em cada dimensão: [xmin, xmax]
    double xl[3] = { 0.0, 0.0, 0.0 };
    double xu[3] = { 1.0, 1.0, 1.0 };

    size_t calls = 500000; // Número de pontos gerados

    gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(3);
    gsl_monte_vegas_integrate(&G, xl, xu, 3, calls, r, s, &result, &error);
    gsl_monte_vegas_free(s);

    std::cout << "Resultado: " << result << " +/- " << error << std::endl;

    gsl_rng_free(r);
    
    }
    */
    
    
                 /* --------------------------- ROOT-------------------------------*/
/*    
    const int fdim = 3; // Número de variáveis (x, y, z)
    int pdim = 3;
    double parr[pdim] = {0.0, 0.0, 0.0}; // Parâmetros adicionais (caso a função dependa de parâmetros)

    double xa[fdim] = {0.0, 0.0, 0.0}; // Limites inferiores da integração
    double xb[fdim] = {1.0, 1.0, 1.0}; // Limites superiores da integração

    double result1 = GetIntegralByRoot(*my_function_ROOT, fdim, parr, pdim, xa, xb,0);

    std::cout << "Resultado: " << result1 << std::endl;
 */
 
    
/*    
    
    if (false){

		int fdim = 2;
		int pdim = 2;
		//                 x   y 
		double a[fdim] = { 0. , 0.};
		double b[fdim] = { 9.2, 3.};

		double parametros[pdim] = {0.,0.};
		//double resultado = GetIntegralByRoot(*integg,fdim,parametros,pdim,a,b); //default 
		double resultado2 = GetIntegralByRoot(*integg,fdim,parametros,pdim,a,b,0); //vegas ROOT
		//double resultado =  GetIntegralMultiDimVegasGSL(a,b,fdim,parametros,pdim,*integg_vegas_gsl); //vegas 
		std::cout << "Resultado gsl vegas = \t"   << resultado2 << "\n";
		
	
	}
    
 */   
    
    
    
                  /*------------------------------------------------------------------------------- */
                 /* --------------------------- Seção de choque ROOT-------------------------------*/
                /*--------------------------------------------------------------------------------*/
                
    
    const int fdim = 2; // Número de variáveis (y3, y4)
    int pdim = 2; // Número de parâmetros

    double xa[fdim] = {0.0, 0.0}; // Limites inferiores da integração
    double xb[fdim] = {1.0, 1.0}; // Limites superiores da integração

    // Abre um arquivo para escrever os dados da integral
    std::ofstream outFile("integral_vs_pT.txt");

    
    // Cria um gráfico TGraph para armazenar os resultados da integral
    TGraph *graph = new TGraph();
    

    
    // Loop para variar pT de 0 a 100
    for (double pT = 0.0; pT <= 10; pT += 0.1) {
        double params[pdim] = {4.5, pT}; // Atualiza o valor de pT nos parâmetros  
        
        //----bottom: 4.5 GeV ---- charm: 1.5 GeV //

        // Calcula a integral para o valor atual de pT
        double resultR = GetIntegralByRoot(*differential_cross_section, fdim, params, pdim, xa, xb, 1); 
        
        graph->SetPoint(graph->GetN(), pT, resultR);
       
        
       // Escreve os resultados no arquivo
        //outFile << pT << "\t" << resultR << std::endl;

        //std::cout << "pT = " << pT << ": Resultado da integral = " << resultR << std::endl;
    }




        // Cria um canvas para plotar o gráfico
        TCanvas *c1 = new TCanvas("c1", "Integral vs  #p_{T}", 900, 500);
        c1->SetGrid();
        
        


        // Define os títulos dos eixos
        graph->SetTitle("Bottom Cross Section");
        graph->GetXaxis()->SetTitle("p_{T} [GeV]");
        graph->GetYaxis()->SetTitle("#frac{#partial#sigma }{#partialp_{T}^{2}}  [mb/GeV^{2}]");

       // Define o tamanho e o estilo dos pontos
        graph->SetMarkerSize(0.5); // Tamanho dos pontos
        graph->SetMarkerStyle(21); // Estilo dos pontos (círculos)
        
       // Cria uma legenda
        TLegend *legend = new TLegend(0.85, 0.85, 0.75, 0.75); // Coordenadas da legenda (ajuste conforme necessário)
        legend->AddEntry(graph, "Teorico", "lp"); // Adiciona a entrada "Dados" com quadrados e linha
        
 
       
       // Plota o gráfico
        graph->Draw("APL");
        legend->Draw();

       // Salva o gráfico em um arquivo (opcional)
        c1->SaveAs("integral_vs_pT.png");


   
    
    
   
    
    
                  /*------------------------------------------------------------------------------- */
                 /* --------------------------- Seção de choque GSL--------------------------------*/
                /*--------------------------------------------------------------------------------*/
    
 /*    

    gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(3); // 3 dimensões

    // Define a função integrand com os parâmetros
    gsl_monte_function integrand = {&differential_cross_section_GSL, 3, 0};

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);

    double xl[3] = {0.0, 0.0, 0.0}; // Limites inferiores (Y3, Y4, PT^2)
    double xu[3] = {1.0, 1.0, 1.0};   // Limites superiores (Y3, Y4, PT^2)

    double resultG, error;
    gsl_monte_vegas_integrate(&integrand, xl, xu, 3, 1000000, rng, state, &resultG, &error);

    std::cout << "Cross Section GSL: " << resultG << " +/- "<< error <<"\n" << "Cross section ROOT:  " << resultR<< std::endl;

    gsl_monte_vegas_free(state);
    gsl_rng_free(rng); */
    
    return 0;
}

