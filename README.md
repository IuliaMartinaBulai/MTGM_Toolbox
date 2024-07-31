# MTGM_Toolbox
MTGM: Metastatic Tumor Growth Modeling Toolbox 

(c) Developed by Iulia Martina Bulai, Maria Carmela De Bonis, Concetta Laurita

In this toolbox we propose an algorithm for computing the biological observables, such as the total metastatic mass (M) and the cumulative number of metastases (N), respectively, from a generalized metastatic tumor growth model without and with treatment, respectively. 
The model describes the primary tumor growth by means of an Ordinary Differential Equation (ODE) and  the evolution of the metastatic density using a transport Partial Differential Equation (PDE), [1]. The numerical computation of the biological observables is based on a Volterra Integral Equation solving method. 

The software implements a numerical method developed in the following very recent paper:
[1] I.M. Bulai, M.C. De Bonis, C. Laurita. Numerical solution of metastatic tumor growth models with treatment. Appl. Math. Comput. (2024) To appear.
Differently than in [1] here we have implemented an iterative method in order to be able to compute the biological observable, M and N respectively, for long interval of times. The toolbox, beside considering the 2D non-autonomous case (metastatic tumor growth with treatment)- kind = 3, introduces also other 4 cases, e.g.
- the case when a simple VIE (Volterra Integral Equation) of the second kind must to be solved, kind = 0
- the 1D model for the metastatic tumor growth model assuming that the solutions of the ODE systems describing the metastatic and the primary tumor growths are known, kind = 1 
- the 1D model for the metastatic tumor growth model assuming that the solutions of the ODE systems describing the metastatic and the primary tumor growths are NOT known, the numerical solutions of such systems will be used, kind = 2 
- the 2D autonomous metastatic tumor growth model, kind = 4   

The algorithms implemented here are described in detail in: 
[2] I.M. Bulai, M.C. De Bonis, C. Laurita. A new MATLAB software for numerical computation of biological observables for metastatic tumor growth. Submitted (2024).


GETTING STARTED
Run

>> mtgm_demo0: This demo computes the weighted solution of a VIE.

>> mtgm_demo1 This demo computes the metastatic mass M(t) and the cumulative number of metastases N(t) for breast tumor data assuming that both the primary and secondary tumors growth laws are Gompertz laws. Moreover for this 1D model the analytical solution of the growth laws is known.

>> mtgm_demo2 Same as mtgm_demo1 with the difference that here the analytical solution of the ODE equation is assumed to be unknown thus it is solved numerically. The results can be compared with those from mtgm_demo1.

>> mtgm_demo3 This demo computes the cumulative number of metastases N(t) and the total metastatic mass M(t) for a 2D non-autonomous metastatic tumor growth model without considering treatment. The relative error corresponding to this simulation is also computed.

>> mtgm_demo4 Same as for mtgm_demo3 but assuming also that three different kinds of therapy are administered. For the sake of brevity no relative error is computed here.

>> mtgm_demo5 This demo computes the cumulative number of metastases N(t) and the total metastatic mass M(t) for a 2D autonomous metastatic tumor growth model. The results can be compared with those from mtgm_demo3.

>> mtgm_demo6 This demo computes the cumulative number of metastases N(t) and the total metastatic mass M(t) for a 2D non-autonomous metastatic tumor growth model considering only one type of treatment and dose but administered with different schedules.

>> mtgm_demo7 This demo computes the cumulative number of metastases N(t) and the total metastatic mass M(t) for a 2D non-autonomous metastatic tumor growth model considering only one type of treatment, same administration schedule but different doses.

>> mtgm_demo8 This demo computes the cumulative number of metastases N(t) and the total metastatic mass M(t) for a 2D non-autonomous metastatic tumor growth model combining two different treatments.

>> mtgm_demo9 Same as for mtgm_demo4 but for T = 360 days instead of T = 15 days.

>> mtgm_demo10 Same as for mtgm_demo3 but for T = 360 days instead of T = 40 days.


-The files argselectAssign.m argselectCheck.m are part of Spectral Graph Wavelet Transform (SGWT)  toolbox and can be downloaded at page: https://wiki.epfl.ch/sgwt
-The file gaussq.m is the result of a translation of the fortran procedure http://www.netlib.org/go/gaussq.f
-The file class1.m . 

License : 

The MTGM toolbox is a MATLAB library released under the GPL.

The MTGM toolbox is free software: you can redistribute it and/or modify it under the terms of the GNU  General Public License as published by the Free Software Foundation, either version 3 of the License,  or (at your option) any later version.

The MTGM toolbox is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the MTTM toolbox. If not, see <http://www.gnu.org/licenses/>.
