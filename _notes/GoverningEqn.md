---
title: "Governing equations of combustion flow"
collection: notes
type: "Note"
permalink: /notes/GoverningEqn
date: 2022-07-19
location: "Beijing, China"
tags: 
   - reactive flow
   - CFD
   - Favre average
---

对带有化学反应流动的控制方程进行了总结，并给出了曲线坐标系下的流动控制方程，用于CFD离散

# 燃烧控制方程推导

## N-S方程

连续方程：

$$
\frac{\partial\rho}{\partial t}+\frac{\partial\rho u_j}{\partial x_j}=0
$$

动量方程(i=1,2,3)：

$$
\frac{\partial\rho u_i}{\partial t}+\frac{\partial\rho u_ju_i}{\partial x_j}=-\frac{\partial p}{\partial x_i}+\frac{\partial\tau_{ji}}{\partial x_j}+F_i
$$

组分方程($N$组分，$k=1,2,\cdots,N$)：

$$
\frac{\partial\rho Y_k}{\partial t}+\frac{\partial\rho u_jY_k}{\partial x_j}=-\frac{\partial \jmath_i^k}{\partial x_j}+\dot{\omega}_k
$$

总焓方程($h_t=h+u_iu_i/2$)：

$$
\frac{\partial\rho h_t}{\partial t}+\frac{\partial\rho u_jh_t}{\partial x_j}=\frac{\partial p}{\partial t}+\frac{\partial (\jmath_j^h+u_i\tau_{ij})}{\partial x_j}+u_jF_j
$$

上述方程中的本构关系：

假设牛顿流体，则有粘性应力张量$\tau_{ij}$：

$$
\tau_{ij}=\mu_l(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i})-\frac23\mu_l\frac{\partial u_k}{\partial x_k}\delta_{ij}
$$

组分分子扩散采用Fick定律，忽略了Soret效应(由于温度梯度导致的组分扩散)以及由于压力梯度导致的组分扩散：

$$
\jmath_j^k=-\rho D_k\frac{\partial Y_k}{\partial x_j}=-\frac{\mu_l}{Sc}\frac{\partial Y_k}{\partial x_j}
$$

焓扩散，忽略了辐射传热和Dufour效应(由于质量分数梯度导致的焓扩散)：

$$
\jmath_j^h=-\frac{\mu_l}{Pr}[\frac{\partial h}{\partial x_j}+\sum_{k=1}^N(\frac{Pr}{Sc_k}-1)h_k\frac{\partial Y_k}{\partial x_j}]
$$

## RANS

### Favre平均

为RANS，由于燃烧中密度通常会发生变化，因此采用雷诺平均会导致密度脉动产生额外的未封闭项，故采用Favre平均。其定义及性质简单推导如下：

$$
\tilde{Q}=\frac{\overline{\rho Q}}{\bar{\rho}}
$$

$$
\bar{\tilde{Q}}=\overline{(\frac{\overline{\rho Q}}{\bar{\rho}})}=\tilde{Q}
$$

$$
\overline{\rho'\tilde{Q}}=\frac1T\int_0^T\rho'\tilde{Q}dt=\frac1T\int_0^T\rho'\frac{\overline{\rho Q}}{\bar{\rho}}dt=\frac{\overline{\rho Q}}{\bar{\rho}}\frac1T\int_0^T\rho'dt=\tilde{Q}\bar{\rho'}=0
$$

$$
\tilde{\tilde{Q}}=\frac{\overline{\rho\tilde{Q}}}{\bar{\rho}}=\frac{\overline{\bar{\rho}\tilde{Q}}+\overline{\rho'\tilde{Q}}}{\bar{\rho}}=\tilde{Q}+\frac{\overline{\rho'\tilde{Q}}}{\bar{\rho}}=\tilde{Q}
$$

### N-S方程的Favre平均

总体的思路都是一样的，对每隔方程取时均，并利用关系$\overline{\rho Q}=\bar{\rho}\tilde{Q}$进行变换，即可得到相应的Favre平均的N-S方程：

连续方程：

$$
\frac{\partial\bar{\rho}}{\partial t}+\frac{\partial{\bar\rho} \tilde{u_j}}{\partial x_j}=0
$$

动量方程(i=1,2,3)：

$$
\overline{\frac{\partial\rho u_i}{\partial t}+\frac{\partial\rho u_ju_i}{\partial x_j}}=\overline{-\frac{\partial p}{\partial x_i}+\frac{\partial\tau_{ji}}{\partial x_j}+F_i}
$$

$$
\frac{\partial\bar{\rho} \tilde{u_i}}{\partial t}+\frac{\partial\bar{\rho} \widetilde{u_ju_i}}{\partial x_j}=-\frac{\partial \bar{p}}{\partial x_i}+\frac{\partial\bar{\tau}_{ji}}{\partial x_j}+\bar{F_i}
$$

**二阶的Favre平均可以化为以下两项**：

$$
\widetilde{u_ju_i}=\widetilde{(\tilde{u}_i+u_i^{''})(\tilde{u}_j+u_j^{''})}=\widetilde{\tilde{u}_i\tilde{u}_j+\tilde{u}_iu_j^{''}+\tilde{u}_ju_i^{''}+u_i^{''}u_j^{''}}=\tilde{u}_i\tilde{u}_j+\widetilde{u_i^{''}u_j^{''}}
$$

因此动量方程可化为：

$$
\frac{\partial\bar{\rho} \tilde{u}_i}{\partial t}+\frac{\partial\bar{\rho}\tilde{u}_i\tilde{u}_j}{\partial x_j}=-\frac{\partial\bar{\rho}\widetilde{u_i^{''}u_j^{''}}}{\partial x_j}-\frac{\partial \bar{p}}{\partial x_i}+\frac{\partial\bar{\tau}_{ji}}{\partial x_j}+\bar{F_i}
$$

同理，组分方程($N$组分，$k=1,2,\cdots,N$)：

$$
\frac{\partial\bar{\rho} \tilde{Y}_k}{\partial t}+\frac{\partial\bar{\rho}\tilde{u}_j\tilde{Y}_k}{\partial x_j}=-\frac{\partial\bar{\rho}\widetilde{u_j^{''}Y_k^{''}}}{\partial x_j}-\frac{\partial \overline{\jmath_i^k}}{\partial x_i}+\bar{\dot{\omega}}_k
$$

总焓方程($h_t=h+u_iu_i/2$)：

$$
\frac{\partial\bar{\rho} \tilde{h}_t}{\partial t}+\frac{\partial\bar{\rho}\tilde{u}_j\tilde{h}_t}{\partial x_j}=-\frac{\partial\bar{\rho}\widetilde{u_j^{''}h_t^{''}}}{\partial x_j}+\frac{\partial\bar{p}}{\partial t}+\frac{\partial(\overline{\jmath_j^h}+\overline{u_i\tau_{ij}})}{\partial x_j}+\overline{u_jF_j}
$$

简单来看，上述方程中的未封闭项有：

1. **雷诺应力$\widetilde{u_i^{''}u_j^{''}}$**：通常直接由非反应流中的湍流模型直接封闭，一般忽略放热对于雷诺应力的影响
2. **组分$\widetilde{u_j^{''}Y_k^{''}}$和温度$\widetilde{u_j^{''}T^{''}}$湍流通量**：通常采用梯度扩散假设进行封闭：
   $$
   \bar{\rho}\widetilde{u_j^{''}Y_k^{''}}=-\frac{\mu_t}{Sc_{kt}}\frac{\partial\tilde{Y_k}}{\partial x_j}
   $$
   但在实际中会出现逆梯度输运
3. **层流扩散通量$\overline{\jmath_i^k}$和$\overline{\jmath_j^h}$**：高雷诺数下，通常认为这都非常小，远小于湍流输运
4. **组分反应率$\bar{\dot{\omega}}_k$**：湍流燃烧模型的目标

# CFD中控制方程推导

## 基本控制方程

直角坐标系下的N-S方程组为
$$
\frac{\partial\boldsymbol{Q}}{\partial t}+\frac{\partial\boldsymbol{f}_c}{\partial x}+\frac{\partial\boldsymbol{g}_c}{\partial y}+\frac{\partial\boldsymbol{h}_c}{\partial z}=\frac{\partial\boldsymbol{f}_v}{\partial x}+\frac{\partial\boldsymbol{g}_v}{\partial y}+\frac{\partial\boldsymbol{h}_v}{\partial z}+\boldsymbol{S}
$$
式中：

$$
\boldsymbol{Q}=\left[\begin{array}{c}
\rho \\\rho u\\\rho v\\\rho w\\\rho E\end{array}
\right],\qquad E=e+\frac12(u^2+v^2+w^2)
$$

$$
\boldsymbol{f}_c=\left[\begin{array}{c}
\rho u\\\rho u^2+p\\\rho uv\\\rho uw\\(\rho E+p)u
\end{array}\right],\quad
\boldsymbol{g}_c=\left[\begin{array}{c}
\rho v\\\rho v^2\\\rho v^2+p\\\rho vw\\(\rho E+p)v
\end{array}\right],\quad
\boldsymbol{h}_c=\left[\begin{array}{c}
\rho w\\\rho uw\\\rho vw\\\rho w^2+p\\(\rho E+p)w
\end{array}\right]
$$

$$
\boldsymbol{f}_v=\left[\begin{array}{c}
0\\\tau_{xx}\\\tau_{xy}\\\tau_{xz}\\b_x
\end{array}\right],\quad
\boldsymbol{g}_v=\left[\begin{array}{c}
0\\\tau_{xy}\\\tau_{yy}\\\tau_{yz}\\b_y
\end{array}\right],\quad
\boldsymbol{h}_v=\left[\begin{array}{c}
0\\\tau_{xz}\\\tau_{yz}\\\tau_{zz}\\b_z
\end{array}\right]
$$

$$
\tau_{ij}=2\mu(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i})-\frac23\mu\frac{\partial u_m}{\partial x_m}\delta_{ij}
$$

$$
b_i=u_j\tau_{ij}-\dot{q}_i
$$

## 无量纲化方程

取下列参考量：
$$
L_{ref}=L_{flow/grid...}\quad T_{ref}=T_\infin \quad\mu_{ref}=\mu_\infin\quad\rho_{ref}=\rho_\infin \quad v_{ref}=V_\infin
$$

$$
p_{ref}=\rho_{ref} v_{ref}^2\quad t_{ref}=L_{ref}/v_{ref}\quad E_{ref}=v_{ref}^2\quad R_{ref}={C_p}_{ref}=\frac{v_{ref}^2}{T_{ref}}
$$

因此式中各项的无量纲表达为
$$
x^*=\frac x{L_{ref}}\quad y^*=\frac y{L_{ref}}\quad z^*=\frac z{L_{ref}}\quad T^*=\frac T{T_{ref}}\quad\mu^*=\frac\mu{\mu_{ref}}\quad\rho^*=\frac\rho{\rho_{ref}}
$$

$$
u^*=\frac u{v_{ref}}\quad v^*=\frac v{v_{ref}}\quad w^*=\frac w{v_{ref}}\quad c^*=\frac c{c_{ref}}
$$

$$
p^*=\frac p{\rho_{ref}v_{ref}^2}\quad t^*=\frac t{L_{ref}/v_{ref}}\quad e^*=\frac e{v_{ref}^2}\quad C_p^*=\frac{C_p}{v_{ref}^2/T_{ref}}\quad R^*=\frac R{v_{ref}^2/T_{ref}}
$$

将变量无量纲化后，得到的无量纲形式控制方程为
$$
\frac{\partial\boldsymbol{Q}}{\partial t}+\frac{\partial\boldsymbol{f}_c}{\partial x}+\frac{\partial\boldsymbol{g}_c}{\partial y}+\frac{\partial\boldsymbol{h}_c}{\partial z}=\frac1{Re_{ref}}(\frac{\partial\boldsymbol{f}_v}{\partial x}+\frac{\partial\boldsymbol{g}_v}{\partial y}+\frac{\partial\boldsymbol{h}_v}{\partial z})+\boldsymbol{S}
$$
式中各项的表达形式相比有量纲形式毫无改变，只是都替换为了相应的无量纲变量，差别在于无量纲方程在扩散项前有参考雷诺数
$$
Re_{ref}=\frac{\rho_\infin V_\infin L_{ref}}{\mu_\infin}
$$
气体状态方程形式也未发生改变：
$$
p^*=\rho^*R^*T^*
$$

## 曲线坐标系形式

方程形式为
$$
\frac{\partial\hat{\boldsymbol{Q}}}{\partial\tau}+\frac{\partial\hat{\boldsymbol{F}}_c}{\partial\xi}+\frac{\partial\hat{\boldsymbol{G}}_c}{\partial\eta}+\frac{\partial\hat{\boldsymbol{H}}_c}{\partial\zeta}=\frac1{Re_{ref}}(\frac{\partial\hat{\boldsymbol{F}}_v}{\partial\xi}+\frac{\partial\hat{\boldsymbol{G}}_v}{\partial\eta}+\frac{\partial\hat{\boldsymbol{H}}_v}{\partial\zeta})+\hat{\boldsymbol{S}}
$$
式中各项与原变量关系如下：
$$
\hat{\boldsymbol{Q}}=\boldsymbol{JQ}\quad \hat{\boldsymbol{S}}=\boldsymbol{JS}
$$
逆变速度$U,V,W$定义如下：

$$
\begin{cases}
   \ U=\xi_xu+\xi_yv+\xi_zw+\xi_t\\
   \ V=\eta_xu+\eta_yv+\eta_zw+\eta_t\\
   \ W=\zeta_xu+\zeta_yv+\zeta_zw+\zeta_t
\end{cases}
$$

通量项定义如下：
$$
\hat{\boldsymbol{F}}_c=\boldsymbol{J}(\xi_x\boldsymbol{f}_c+\xi_y\boldsymbol{g}_c+\xi_z\boldsymbol{h}_c+\xi_t\boldsymbol{Q})=\boldsymbol{J}\left[\begin{array}{c}\rho U\\\rho uU+p\xi_x\\\rho vU+p\xi_y\\\rho wU+p\xi_z\\(\rho E+p)U-p\xi_t
\end{array}\right]
$$

$$
\hat{\boldsymbol{F}}_v=\boldsymbol{J}(\xi_x\boldsymbol{f}_v+\xi_y\boldsymbol{g}_v+\xi_z\boldsymbol{h}_v)=\boldsymbol{J}\left[\begin{array}{c}0\\\xi_x\tau_{xx}+\xi_y\tau_{xy}+\xi_z\tau_{xz}\\\xi_x\tau_{xy}+\xi_y\tau_{yy}+\xi_z\tau_{yz}\\\xi_x\tau_{xz}+\xi_y\tau_{yz}+\xi_z\tau_{zz}\\\xi_xb_x+\xi_yb_y+\xi_zb_z
\end{array}\right]
$$

$$
\hat{\boldsymbol{G}}_c=\boldsymbol{J}(\eta_x\boldsymbol{f}_c+\eta_y\boldsymbol{g}_c+\eta_z\boldsymbol{h}_c+\eta_t\boldsymbol{Q})=\boldsymbol{J}\left[\begin{array}{c}\rho V\\\rho uV+p\eta_x\\\rho vV+p\eta_y\\\rho wV+p\eta_z\\(\rho E+p)V-p\eta_t
\end{array}\right]
$$

$$
\hat{\boldsymbol{G}}_v=\boldsymbol{J}(\eta_x\boldsymbol{f}_v+\eta_y\boldsymbol{g}_v+\eta_z\boldsymbol{h}_v)=\boldsymbol{J}\left[\begin{array}{c}0\\\eta_x\tau_{xx}+\eta_y\tau_{xy}+\eta_z\tau_{xz}\\\eta_x\tau_{xy}+\eta_y\tau_{yy}+\eta_z\tau_{yz}\\\eta_x\tau_{xz}+\eta_y\tau_{yz}+\eta_z\tau_{zz}\\\eta_xb_x+\eta_yb_y+\eta_zb_z
\end{array}\right]
$$

$$
\hat{\boldsymbol{H}}_c=\boldsymbol{J}(\zeta_x\boldsymbol{f}_c+\zeta_y\boldsymbol{g}_c+\zeta_z\boldsymbol{h}_c+\zeta_t\boldsymbol{Q})=\boldsymbol{J}\left[\begin{array}{c}\rho W\\\rho uW+p\zeta_x\\\rho vW+p\zeta_y\\\rho wW+p\zeta_z\\(\rho E+p)W-p\zeta_t
\end{array}\right]
$$

$$
\hat{\boldsymbol{H}}_v=\boldsymbol{J}(\zeta_x\boldsymbol{f}_v+\zeta_y\boldsymbol{g}_v+\zeta_z\boldsymbol{h}_v)=\boldsymbol{J}\left[\begin{array}{c}0\\\zeta_x\tau_{xx}+\zeta_y\tau_{xy}+\zeta_z\tau_{xz}\\\zeta_x\tau_{xy}+\zeta_y\tau_{yy}+\zeta_z\tau_{yz}\\\zeta_x\tau_{xz}+\zeta_y\tau_{yz}+\zeta_z\tau_{zz}\\\zeta_xb_x+\zeta_yb_y+\zeta_zb_z
\end{array}\right]
$$

以上即为曲线坐标形式的N-S方程
