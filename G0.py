 
import streamlit as st
from iapws import IAPWS97
import numpy as np
import pandas as pd
import math as M
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
 

def iso_bar(wsp_point, min_s=-0.1, max_s=0.11, step_s=0.011, color='r'):
    if not isinstance(wsp_point, list):
        iso_bar_0_s = np.arange(wsp_point.s + min_s, wsp_point.s + max_s, step_s).tolist()
        iso_bar_0_h = [WSP(P=wsp_point.P, s=i).h for i in iso_bar_0_s]
    else:
        iso_bar_0_s = np.arange(wsp_point[0].s + min_s, wsp_point[1].s + max_s, step_s).tolist()
        iso_bar_0_h = [WSP(P=wsp_point[1].P, s=i).h for i in iso_bar_0_s]
    plt.plot(iso_bar_0_s, iso_bar_0_h, color)


st.write('# Курсовая работа Part II')

st.sidebar.header("Ввод параметров:")

st.sidebar.markdown('**1) Средний диаметр:**')
value_d = st.sidebar.slider('',60, 150, (90, 111))
st.sidebar.write('Значения cредних диаметров: ', value_d[0], "-", value_d[1], "м")

st.sidebar.markdown('**2) Давление пара перед турбиной:**')
P_0 = st.sidebar.number_input('', value=29)
st.sidebar.write('P_0 = ', P_0, "МПа")

st.sidebar.markdown('**3) Температура пара перед турбиной**')
t_0 = st.sidebar.number_input('', value=575)
st.sidebar.write('t_0 = ', t_0, "°C")

st.sidebar.markdown('**4) Частота вращения ротора турбины**')
n_ = st.sidebar.number_input('', value=50)
st.sidebar.write('n = ',n_, "c^(-1)")


st.sidebar.markdown('**5) Расход водяного пара**')
G = st.sidebar.number_input('', value=448.34)
st.sidebar.write('G_0 = ', G, "кг/с")


st.sidebar.markdown('**6) Располагаемый теплоперепад ступени**')
H = st.sidebar.number_input('', value=110)
st.sidebar.write('H_0 = ', H, "кДж/кг")

def iso_bar(wsp_point, min_s=-0.1, max_s=0.11, step_s=0.011, color = 'b'):
    if not isinstance(wsp_point,list):
        iso_bar_0_s = np.arange(wsp_point.s+min_s,wsp_point.s+max_s,step_s).tolist()
        iso_bar_0_h = [IAPWS97(P = wsp_point.P, s = i).h for i in iso_bar_0_s]
    else:
        iso_bar_0_s = np.arange(wsp_point[0].s+min_s,wsp_point[1].s+max_s,step_s).tolist()
        iso_bar_0_h = [IAPWS97(P = wsp_point[1].P, s = i).h for i in iso_bar_0_s]
    plt.plot(iso_bar_0_s,iso_bar_0_h,color)
    
d = 1.1 #m
p_0 = 28 #МПа
t_0 = 570 #град Цельсия
T_0 = t_0+273.15 #K
n = 60 #Гц
G_0 = 443.936265 #кг/с
G_k = 292.733965 #кг/с
H_0 = 100 #кДж/кг
rho = 0.05 #степень реактивности
l_1 = 0.015 #м
alpha_1 = 12 #град
b_1 = 0.06 #м
Delta = 0.003 #м
b_2 = 0.03 #м
kappa_vs = 0 #коэффициент использования выходной скорости


D1 = drs - deltaD
sat_steam = IAPWS97(P=P0, h=h0)
s_0 = sat_steam.s
t_0 = sat_steam.T
error = 2
i = 1
while error > 0.5:
    rho = rho_s + 1.8 / (tetta + 1.8)
    X = (fi * M.cos(M.radians(alfa))) / (2 * M.sqrt(1 - rho))
    H01 = 12.3 * (D1 / X) ** 2 * (n / 50) ** 2
    h2t = h0 - H01
    steam2t = IAPWS97(h=h2t, s=s_0)
    v2t = steam2t.v
    l11 = G0 * v2t * X / (M.pi ** 2 * D1 ** 2 * n * M.sqrt(1 - rho) * M.sin(M.radians(alfa)) * mu1)
    tetta_old = tetta
    tetta = D1 / l11
    error = abs(tetta - tetta_old) / tetta_old * 100
    i += 1

l21 = l11 + delta
d_s = D1 - l21
steam_tz = IAPWS97(P=Pz, s=s_0)
h_zt = steam_tz.h
H0 = h0 - h_zt
Hi = H0 * etaoi
h_z = h0 - Hi
steam_z = IAPWS97(P=Pz, h=h_z)
v_2z = steam_z.v
x = Symbol('x')
с = solve(x ** 2 + x * d_s - (l21 * (d_s + l21) * v_2z / v2t))
for j in с:
    if j > 0:
        l2z = j
d2z = d_s + l2z
tetta1 = (l21 + d_s) / l21
tettaz = (l2z + d_s) / l2z
rho1 = rho_s + 1.8 / (1.8 + tetta1)
rhoz = rho_s + 1.8 / (1.8 + tettaz)
X1 = (fi * M.cos(M.radians(alfa))) / (2 * M.sqrt(1 - rho1))
Xz = (fi * M.cos(M.radians(alfa))) / (2 * M.sqrt(1 - rhoz))

DeltaZ = 1
ite = 0
while DeltaZ > 0:
    matr = []
    Num = 0
    SumH = 0
    for _ in range(int(Z)):
        li = (l21 - l2z) / (1 - Z) * Num + l21
        di = (D1 - d2z) / (1 - Z) * Num + D1
        tetta_i = di / li
        rho_i = rho_s + 1.8 / (1.8 + tetta_i)
        X_i = (fi * M.cos(M.radians(alfa))) / (2 * M.sqrt(1 - rho_i))
        if Num < 1:
            H_i = 12.3 * (di / X_i) ** 2 * (n / 50) ** 2
        else:
            H_i = 12.3 * (di / X_i) ** 2 * (n / 50) ** 2 * 0.95
        Num = Num + 1
        H_d = 0
        SumH = SumH + H_i
        matr.append([Num, round(di, 3), round(li, 3), round(tetta_i, 2), round(rho_i, 3), round(X_i, 3), round(H_i, 2),
                     round(H_d, 2)])
    H_m = SumH / Z
    q_t = 4.8 * 10 ** (-4) * (1 - etaoi) * H0 * (Z - 1) / Z
    Z_new = round(H0 * (1 + q_t) / H_m)
    DeltaZ = abs(Z - Z_new)
    # print(ite, Z)
    Z = Z_new
    ite += 1
DeltaH = (H0 * (1 + q_t) - SumH) / Z
a = 0
for elem in matr:
    matr[a][7] = round(elem[6] + DeltaH, 2)
    a += 1

N_ = []
di_ = []
li_ = []
tettai_ = []
rhoi_ = []
Xi_ = []
Hi_ = []
Hdi_ = []
a = 0
for elem in matr:
    N_.append(matr[a][0])
    di_.append(matr[a][1])
    li_.append(matr[a][2])
    tettai_.append(matr[a][3])
    rhoi_.append(matr[a][4])
    Xi_.append(matr[a][5])
    Hi_.append(matr[a][6])
    Hdi_.append(matr[a][7])
    a += 1

di_ = [float(x) for x in di_]
li_ = [float(x) for x in li_]
tettai_ = [float(x) for x in tettai_]
rhoi_ = [float(x) for x in rhoi_]
Xi_ = [float(x) for x in Xi_]
Hi_ = [float(x) for x in Hi_]
Hdi_ = [float(x) for x in Hdi_]

## Таблица
table = pd.DataFrame({"№ ступени": (N_),
                      "di, м": (di_),
                      "li, м": (li_),
                      "θi ": (tettai_),
                      "ρi ": (rhoi_),
                      "Xi ": (Xi_),
                      "Hi, кДж/кг": (Hi_),
                      "Hi + Δ, кДж/кг": (Hdi_)
                      }
                     )

st.dataframe(table)

## Графики
z = []
for a in range(1, Z + 1):
    z.append(a)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 15, 1))
plt.grid(True)
plt.plot(z, di_, '-ro')
plt.title('Рисунок 1. Распределение средних диаметров по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 15, 1))
plt.grid(True)
plt.plot(z, li_, '-ro')
plt.title('Рисунок 2. Распределение высот лопаток по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 15, 1))
plt.grid(True)
plt.plot(z, tettai_, '-ro')
plt.title('Рисунок 3. Распределение обратной веерности по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 15, 1))
plt.grid(True)
plt.plot(z, rhoi_, '-ro')
plt.title('Рисунок 4. Распределение степени реактивности по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 15, 1))
plt.grid(True)
plt.plot(z, Xi_, '-ro')
plt.title('Рисунок 5. Распределение U/Cф по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 15, 1))
plt.grid(True)
plt.plot(z, Hi_, '-ro')
plt.title('Рисунок 6. Распределение теплоперепадов по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 15, 1))
plt.grid(True)
plt.plot(z, Hdi_, '-ro')
plt.title('Рисунок 7. Распределение теплоперепадов с учетом невязки по проточной части')
st.pyplot(fig)





