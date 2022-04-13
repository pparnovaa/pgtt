 
import streamlit as st
from iapws import IAPWS97
import numpy as np
import pandas as pd
import matplotlib as plt
import matplotlib.pyplot as plt
st.write('# Курсовая работа по курсу "Паровые и газовые турбины"')
st.write('# Выполнила: Парнова Екатерина')
Ne = 816e6 #МВт
p0 = 12.5e6 #МПа
t0 = 552 #C
T0 = t0 + 273.15 #K 
ppp = 3.74e6 #МПа
tpp = 558 #C
Tpp = tpp+273.15 #K
pk = list(range(int(2e3),int(10e3),500))
tpv = 275 #C
Tpv = tpv+273.15 #K

pk_min = 2e3
pk_max = 10e3

delta_p_0 = 0.05*p0
delta_p_pp = 0.08*ppp
delta_p = 0.03*ppp

from bokeh.plotting import figure
def Calculate_G0_Gk(N_e, p_0, T_0, p_pp, T_pp, p_k, T_pv):
    # Потери:
    d_p0 = 0.05
    d_p_pp = 0.1
    d_p = 0.03
    # Параметры свежего пара
    point_0 = IAPWS97(P=p_0*10**(-6),T=T_0)
    s_0 = point_0.s
    h_0 = point_0.h
    v_0 = point_0.v
    # 
    p_0_ = p_0-0.05*p_0
    point_p_0_ = IAPWS97(P=p_0_*10**(-6),h=h_0)
    t_0_ = point_p_0_.T-273.15
    s_0_ = point_p_0_.s
    v_0_ = point_p_0_.v
    #Теоретический процесс расширения в ЦВД
    p_1t = p_pp+0.1*p_pp
    point_1t = IAPWS97(P=p_1t*10**(-6),s=s_0)
    t_1t = point_1t.T-273.15
    h_1t = point_1t.h
    v_1t = point_1t.v
    #
    point_pp = IAPWS97(P=p_pp*10**(-6),T=T_pp)
    h_pp = point_pp.h
    s_pp = point_pp.s
    v_pp = point_pp.v
    #Действительный процесс расширения в ЦВД
    H_0 = h_0-h_1t
    eta_oi = 0.85
    H_i_cvd = H_0*eta_oi
    h_1 = h_0 - H_i_cvd
    point_1 = IAPWS97(P = p_1t*10**(-6),h = h_1)
    s_1 = point_1.s
    T_1 = point_1.T
    v_1 = point_1.v
    #
    p_pp_ = p_pp - 0.03*p_pp
    point_pp_ = IAPWS97(P=p_pp_*10**(-6),h = h_pp)
    s_pp_ = point_pp_.s
    v_pp_ = point_pp_.v
    #
    point_kt = IAPWS97(P = p_k*10**(-6),s = s_pp)
    T_kt = point_kt.T
    h_kt = point_kt.h
    v_kt = point_kt.v
    s_kt = s_pp
    #
    H_0_csdcnd = h_pp-h_kt
    eta_oi = 0.85
    H_i_csdcnd = H_0_csdcnd*eta_oi
    h_k = h_pp - H_i_csdcnd
    point_k = IAPWS97(P = p_k*10**(-6), h = h_k)
    T_k = point_k.T
    s_k = point_k.s
    v_k = point_k.v
    #
    point_k_v = IAPWS97(P = p_k*10**(-6),x=0)
    h_k_v = point_k_v.h
    s_k_v = point_k_v.s
    eta_oiI = (h_1-h_0)/(h_1t-h_0)
    p_pv = 1.4*p_0
    point_pv = IAPWS97(P = p_pv*10**(-6),T=T_pv)
    h_pv = point_pv.h
    s_pv = point_pv.s
    #
    ksi_pp_oo = 1-(1-(T_k*(s_pp-s_k_v))/((h_0-h_1t)+(h_pp-h_k_v)))/(1-(T_k*(s_pp-s_pv))/((h_0-h_1t)+(h_pp-h_pv)))
    T_0_= IAPWS97(P = p_pv*10**(-6),x = 0).T
    T_ = (point_pv.T - point_k.T) / (T_0_ - point_k.T)
    if T_<=0.636364:
      ksi1=-1.53*T_**2+2.1894*T_+0.0048
    elif 0.636364<T_<=0.736364:
      ksi1=-1.3855*T_**2+2.00774*T_+0.0321
    elif 0.736364<T_<=0.863636:
      ksi1=-2.6536*T_**2+4.2556*T_-0.8569
    
    if T_ <=0.631818 :
        ksi2 = -1.7131*T_**2+2.3617*T_-0.0142
    elif 0.631818<T_<=0.718182:
        ksi2 = -2.5821*T_**2+3.689*T_-0.4825
    elif 0.718182<T_<=0.827273:
        ksi2 =-2.5821*T_**2+3.138*T_-0.3626
          
    ksi = (ksi1+ksi2)/2
    ksi_r_pp = ksi*ksi_pp_oo
    eta_ir = (H_i_cvd+H_i_csdcnd)/(H_i_cvd+(h_pp-h_k_v))*1/(1-ksi_r_pp)
    H_i = eta_ir*((h_0-h_pv)+(h_pp-h_1))
    eta_m = 0.994
    eta_eg = 0.99
    G_0 = N_e/(H_i*eta_m*eta_eg*(10**3))
    G_k = N_e/((h_k-h_k_v)*eta_m*eta_eg*(10**3))*(1/eta_ir-1)
    return eta_ir
eta=[Calculate_G0_Gk(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = p, T_pv = Tpv) for p in pk]

def Calculate_G0(N_e, p_0, T_0, p_pp, T_pp, p_k, T_pv):
    # Потери:
    d_p0 = 0.05
    d_p_pp = 0.1
    d_p = 0.03
    # Параметры свежего пара
    point_0 = IAPWS97(P=p_0*10**(-6),T=T_0)
    s_0 = point_0.s
    h_0 = point_0.h
    v_0 = point_0.v
    # 
    p_0_ = p_0-0.05*p_0
    point_p_0_ = IAPWS97(P=p_0_*10**(-6),h=h_0)
    t_0_ = point_p_0_.T-273.15
    s_0_ = point_p_0_.s
    v_0_ = point_p_0_.v
    #Теоретический процесс расширения в ЦВД
    p_1t = p_pp+0.1*p_pp
    point_1t = IAPWS97(P=p_1t*10**(-6),s=s_0)
    t_1t = point_1t.T-273.15
    h_1t = point_1t.h
    v_1t = point_1t.v
    #
    point_pp = IAPWS97(P=p_pp*10**(-6),T=T_pp)
    h_pp = point_pp.h
    s_pp = point_pp.s
    v_pp = point_pp.v
    #Действительный процесс расширения в ЦВД
    H_0 = h_0-h_1t
    eta_oi = 0.85
    H_i_cvd = H_0*eta_oi
    h_1 = h_0 - H_i_cvd
    point_1 = IAPWS97(P = p_1t*10**(-6),h = h_1)
    s_1 = point_1.s
    T_1 = point_1.T
    v_1 = point_1.v
    #
    p_pp_ = p_pp - 0.03*p_pp
    point_pp_ = IAPWS97(P=p_pp_*10**(-6),h = h_pp)
    s_pp_ = point_pp_.s
    v_pp_ = point_pp_.v
    #
    point_kt = IAPWS97(P = p_k*10**(-6),s = s_pp)
    T_kt = point_kt.T
    h_kt = point_kt.h
    v_kt = point_kt.v
    s_kt = s_pp
    #
    H_0_csdcnd = h_pp-h_kt
    eta_oi = 0.85
    H_i_csdcnd = H_0_csdcnd*eta_oi
    h_k = h_pp - H_i_csdcnd
    point_k = IAPWS97(P = p_k*10**(-6), h = h_k)
    T_k = point_k.T
    s_k = point_k.s
    v_k = point_k.v
    #
    point_k_v = IAPWS97(P = p_k*10**(-6),x=0)
    h_k_v = point_k_v.h
    s_k_v = point_k_v.s
    eta_oiI = (h_1-h_0)/(h_1t-h_0)
    p_pv = 1.4*p_0
    point_pv = IAPWS97(P = p_pv*10**(-6),T=T_pv)
    h_pv = point_pv.h
    s_pv = point_pv.s
    #
    ksi_pp_oo = 1-(1-(T_k*(s_pp-s_k_v))/((h_0-h_1t)+(h_pp-h_k_v)))/(1-(T_k*(s_pp-s_pv))/((h_0-h_1t)+(h_pp-h_pv)))
    T_0_= IAPWS97(P = p_pv*10**(-6),x = 0).T
    T_ = (point_pv.T - point_k.T) / (T_0_ - point_k.T)
    if T_<=0.636364:
      ksi1=-1.53*T_**2+2.1894*T_+0.0048
    elif 0.636364<T_<=0.736364:
      ksi1=-1.3855*T_**2+2.00774*T_+0.0321
    elif 0.736364<T_<=0.863636:
      ksi1=-2.6536*T_**2+4.2556*T_-0.8569
    
    if T_ <=0.631818 :
        ksi2 = -1.7131*T_**2+2.3617*T_-0.0142
    elif 0.631818<T_<=0.718182:
        ksi2 = -2.5821*T_**2+3.689*T_-0.4825
    elif 0.718182<T_<=0.827273:
        ksi2 =-2.5821*T_**2+3.138*T_-0.3626
    ksi = (ksi1+ksi2)/2
    ksi_r_pp = ksi*ksi_pp_oo
    eta_ir = (H_i_cvd+H_i_csdcnd)/(H_i_cvd+(h_pp-h_k_v))*1/(1-ksi_r_pp)
    H_i = eta_ir*((h_0-h_pv)+(h_pp-h_1))
    eta_m = 0.994
    eta_eg = 0.99
    G_0 = N_e/(H_i*eta_m*eta_eg*(10**3))
    G_k = N_e/((h_k-h_k_v)*eta_m*eta_eg*(10**3))*(1/eta_ir-1)
    return G_0
G_0=[Calculate_G0(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = p, T_pv = Tpv) for p in pk]

def Calculate_Gk(N_e, p_0, T_0, p_pp, T_pp, p_k, T_pv):
    # Потери:
    d_p0 = 0.05
    d_p_pp = 0.1
    d_p = 0.03
    # Параметры свежего пара
    point_0 = IAPWS97(P=p_0*10**(-6),T=T_0)
    s_0 = point_0.s
    h_0 = point_0.h
    v_0 = point_0.v
    # 
    p_0_ = p_0-0.05*p_0
    point_p_0_ = IAPWS97(P=p_0_*10**(-6),h=h_0)
    t_0_ = point_p_0_.T-273.15
    s_0_ = point_p_0_.s
    v_0_ = point_p_0_.v
    #Теоретический процесс расширения в ЦВД
    p_1t = p_pp+0.1*p_pp
    point_1t = IAPWS97(P=p_1t*10**(-6),s=s_0)
    t_1t = point_1t.T-273.15
    h_1t = point_1t.h
    v_1t = point_1t.v
    #
    point_pp = IAPWS97(P=p_pp*10**(-6),T=T_pp)
    h_pp = point_pp.h
    s_pp = point_pp.s
    v_pp = point_pp.v
    #Действительный процесс расширения в ЦВД
    H_0 = h_0-h_1t
    eta_oi = 0.85
    H_i_cvd = H_0*eta_oi
    h_1 = h_0 - H_i_cvd
    point_1 = IAPWS97(P = p_1t*10**(-6),h = h_1)
    s_1 = point_1.s
    T_1 = point_1.T
    v_1 = point_1.v
    #
    p_pp_ = p_pp - 0.03*p_pp
    point_pp_ = IAPWS97(P=p_pp_*10**(-6),h = h_pp)
    s_pp_ = point_pp_.s
    v_pp_ = point_pp_.v
    #
    point_kt = IAPWS97(P = p_k*10**(-6),s = s_pp)
    T_kt = point_kt.T
    h_kt = point_kt.h
    v_kt = point_kt.v
    s_kt = s_pp
    #
    H_0_csdcnd = h_pp-h_kt
    eta_oi = 0.85
    H_i_csdcnd = H_0_csdcnd*eta_oi
    h_k = h_pp - H_i_csdcnd
    point_k = IAPWS97(P = p_k*10**(-6), h = h_k)
    T_k = point_k.T
    s_k = point_k.s
    v_k = point_k.v
    #
    point_k_v = IAPWS97(P = p_k*10**(-6),x=0)
    h_k_v = point_k_v.h
    s_k_v = point_k_v.s
    eta_oiI = (h_1-h_0)/(h_1t-h_0)
    p_pv = 1.4*p_0
    point_pv = IAPWS97(P = p_pv*10**(-6),T=T_pv)
    h_pv = point_pv.h
    s_pv = point_pv.s
    #
    ksi_pp_oo = 1-(1-(T_k*(s_pp-s_k_v))/((h_0-h_1t)+(h_pp-h_k_v)))/(1-(T_k*(s_pp-s_pv))/((h_0-h_1t)+(h_pp-h_pv)))
    T_0_= IAPWS97(P = p_pv*10**(-6),x = 0).T
    T_ = (point_pv.T - point_k.T) / (T_0_ - point_k.T)
    if T_<=0.636364:
      ksi1=-1.53*T_**2+2.1894*T_+0.0048
    elif 0.636364<T_<=0.736364:
      ksi1=-1.3855*T_**2+2.00774*T_+0.0321
    elif 0.736364<T_<=0.863636:
      ksi1=-2.6536*T_**2+4.2556*T_-0.8569
    
    if T_ <=0.631818 :
        ksi2 = -1.7131*T_**2+2.3617*T_-0.0142
    elif 0.631818<T_<=0.718182:
        ksi2 = -2.5821*T_**2+3.689*T_-0.4825
    elif 0.718182<T_<=0.827273:
        ksi2 =-2.5821*T_**2+3.138*T_-0.3626
    ksi = (ksi1+ksi2)/2
    ksi_r_pp = ksi*ksi_pp_oo
    eta_ir = (H_i_cvd+H_i_csdcnd)/(H_i_cvd+(h_pp-h_k_v))*1/(1-ksi_r_pp)
    H_i = eta_ir*((h_0-h_pv)+(h_pp-h_1))
    eta_m = 0.994
    eta_eg = 0.99
    G_0 = N_e/(H_i*eta_m*eta_eg*(10**3))
    G_k = N_e/((h_k-h_k_v)*eta_m*eta_eg*(10**3))*(1/eta_ir-1)
    return G_k
Gk=[Calculate_Gk(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = p, T_pv = Tpv) for p in pk]

itog=pd.DataFrame({
"Давление в конденсаторе": (list(range(2000,10000,500))),
"КПД": [Calculate_G0_Gk(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = p, T_pv = Tpv) for p in pk],
"G_0": [Calculate_G0(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = p, T_pv = Tpv) for p in pk],
"G_k": [Calculate_Gk(N_e = Ne, p_0 = p0, T_0 = T0, p_pp = ppp, T_pp = Tpp, p_k = p, T_pv = Tpv) for p in pk]
})


x = (list(range(2000,10000,500)))
y = (eta)

p = figure(
title='Зависимость КПД от давления в конденсаторе',
x_axis_label='давление в конденсаторе',
y_axis_label='КПД')

p.line(x, y, legend_label='Зависимость КПД от давления в конденсаторе', line_width=4)
st.bokeh_chart(p, use_container_width=True)

fighs = plt.figure()
point_0 = IAPWS97(P=p0*1e-6, T=T0)
p_0_d = p0 - delta_p_0
point_0_d = IAPWS97(P=p_0_d*1e-6, h=point_0.h)
p_1t = ppp + delta_p_pp
point_1t = IAPWS97(P=p_1t*10**(-6), s=point_0.s)
H_01 = point_0.h - point_1t.h
kpd_oi = 0.85
H_i_cvd = H_01 * kpd_oi
h_1 = point_0.h - H_i_cvd
point_1 = IAPWS97(P=p_1t*1e-6, h=h_1)
point_pp = IAPWS97(P=ppp*1e-6, T=Tpp)
p_pp_d = ppp - delta_p_pp
point_pp_d = IAPWS97(P=p_pp_d*1e-6, h=point_pp.h)
point_kt = IAPWS97(P=pk_min*1e-6, s=point_pp.s)
H_02 = point_pp.h - point_kt.h
kpd_oi = 0.85
H_i_csd_cnd = H_02 * kpd_oi
h_k = point_pp.h - H_i_csd_cnd
point_k = IAPWS97(P=pk_min*1e-6, h=h_k)

s_0 = [point_0.s-0.05,point_0.s,point_0.s+0.05]
h_0 = [IAPWS97(P = p0*1e-6,s = s_).h for s_ in s_0]
s_1 = [point_0.s-0.05,point_0.s,point_0.s+0.18]
h_1 = [IAPWS97(P=p_1t*1e-6, s = s_).h for s_ in s_1]
s_0_d = [point_0_d.s-0.05, point_0_d.s, point_0_d.s+0.05]
h_0_d = h_0
s_pp = [point_pp.s-0.05,point_pp.s,point_pp.s+0.05]
h_pp = [IAPWS97(P=ppp*1e-6, s=s_).h for s_ in s_pp]
s_k = [point_pp.s-0.05,point_pp.s,point_pp.s+0.8]
h_k = [IAPWS97(P=pk_min*1e-6, s=s_).h for s_ in s_k]
s_pp_d = [point_pp_d.s-0.05,point_pp_d.s,point_pp_d.s+0.05]
h_pp_d = h_pp

plt.plot([point_0.s,point_0.s,point_0_d.s,point_1.s],[point_1t.h,point_0.h,point_0.h,point_1.h],'-or')
plt.plot([point_pp.s,point_pp.s,point_pp_d.s,point_k.s],[point_kt.h,point_pp.h,point_pp.h,point_k.h],'-or')
plt.plot(s_0,h_0)
plt.plot(s_1,h_1)
plt.plot(s_0_d,h_0_d)
plt.plot(s_pp,h_pp)
plt.plot(s_k,h_k)
plt.plot(s_pp_d,h_pp_d)

for x, y, ind in zip([point_pp.s, point_k.s], [point_pp.h, point_k.h], ['{пп}', '{к}']):
  plt.text(x-0.45, y+40, '$h_' + ind + ' = %.2f $'%y)
for x, y, ind in zip([point_kt.s, point_pp_d.s], [point_kt.h, point_pp_d.h], ['{кт}', '{ппд}']):
  plt.text(x+0.03, y+40, '$h_' + ind + ' = %.2f $'%y)

for x, y, ind in zip ([point_0.s, point_1.s], [point_0.h, point_1.h], ['{0}', '{1}']):
  plt.text(x-0.01, y+120, '$h_' + ind + ' = %.2f $'%y)

for x, y, ind in zip([point_1t.s, point_0_d.s], [point_1t.h, point_0_d.h], ['{1т}', '{0д}']):
  plt.text(x+0.03, y-60, '$h_' + ind + ' = %.2f $'%y)


plt.title("h - s диаграмма")
plt.xlabel("Значение энтропии s, кДж/(кг*С)")
plt.ylabel("Значение энтальпии h, кДж/кг")
plt.grid(True)
st.pyplot(fighs)

itog

st.write("Максимальный КПД:")
itog.iloc[0:1]
