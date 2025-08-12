# -*- coding: utf-8 -*-
#s.gonzalezt2
import pulp as lp

#-----------------
#Conjuntos
#-----------------
M=["Luria-Bertani (LB)", "Tripticaseína Soya Agar (TSA)", "Man-Rogosa-Sharpe (MRS)", "Sabouraud Dextrosa Agar (SDA)"]
N=["Agar", "Levadura", "Peptona", "NaCl", "Fosfato"]
R=["fosfato", "peptidos", "aminoacidos", "extracto de levadura", "agar"]
#-----------------
#Parámetros
#-----------------
J={("Luria-Bertani (LB)","Agar"):8.4, ("Tripticaseína Soya Agar (TSA)","Agar"):4.33, ("Man-Rogosa-Sharpe (MRS)","Agar"): 1.2, ("Sabouraud Dextrosa Agar (SDA)","Agar"):0,
   ("Luria-Bertani (LB)","Levadura"):1.26, ("Tripticaseína Soya Agar (TSA)","Levadura"):1.42, ("Man-Rogosa-Sharpe (MRS)","Levadura"): 0, ("Sabouraud Dextrosa Agar (SDA)","Levadura"):0,
   ("Luria-Bertani (LB)","Peptona"):10.22, ("Tripticaseína Soya Agar (TSA)","Peptona"):5.3, ("Man-Rogosa-Sharpe (MRS)","Peptona"): 0.8, ("Sabouraud Dextrosa Agar (SDA)","Peptona"):0.5,
   ("Luria-Bertani (LB)","NaCl"):3, ("Tripticaseína Soya Agar (TSA)","NaCl"):1.26, ("Man-Rogosa-Sharpe (MRS)","NaCl"): 0.025, ("Sabouraud Dextrosa Agar (SDA)","NaCl"):0,
   ("Luria-Bertani (LB)","Fosfato"):0, ("Tripticaseína Soya Agar (TSA)","Fosfato"):0, ("Man-Rogosa-Sharpe (MRS)","Fosfato"): 0.8, ("Sabouraud Dextrosa Agar (SDA)","Fosfato"):0}


G={("Luria-Bertani (LB)","Agar"):19, ("Tripticaseína Soya Agar (TSA)","Agar"):26.1, ("Man-Rogosa-Sharpe (MRS)","Agar"): 11, ("Sabouraud Dextrosa Agar (SDA)","Agar"):25,
   ("Luria-Bertani (LB)","Levadura"):22.8, ("Tripticaseína Soya Agar (TSA)","Levadura"):15.4, ("Man-Rogosa-Sharpe (MRS)","Levadura"): 0, ("Sabouraud Dextrosa Agar (SDA)","Levadura"):0,
   ("Luria-Bertani (LB)","Peptona"):36, ("Tripticaseína Soya Agar (TSA)","Peptona"):46, ("Man-Rogosa-Sharpe (MRS)","Peptona"): 11.4, ("Sabouraud Dextrosa Agar (SDA)","Peptona"):74,
   ("Luria-Bertani (LB)","NaCl"):10, ("Tripticaseína Soya Agar (TSA)","NaCl"):10, ("Man-Rogosa-Sharpe (MRS)","NaCl"): 1, ("Sabouraud Dextrosa Agar (SDA)","NaCl"):0,
   ("Luria-Bertani (LB)","Fosfato"):0, ("Tripticaseína Soya Agar (TSA)","Fosfato"):0, ("Man-Rogosa-Sharpe (MRS)","Fosfato"): 72, ("Sabouraud Dextrosa Agar (SDA)","Fosfato"):0}

H={("Agar", "fosfato"):0, ("Levadura","fosfato"):0, ("Peptona","fosfato"):0, ("NaCl","fosfato"):0, ("Fosfato","fosfato"):100,
   ("Agar", "peptidos"):0, ("Levadura","peptidos"):0, ("Peptona","peptidos"):70, ("NaCl","peptidos"):0, ("Fosfato","peptidos"):0,
   ("Agar", "aminoacidos"):0, ("Levadura","aminoacidos"):0, ("Peptona","aminoacidos"):48, ("NaCl","aminoacidos"):0, ("Fosfato","aminoacidos"):52,
   ("Agar", "extracto de levadura"):0, ("Levadura","extracto de levadura"):76, ("Peptona","extracto de levadura"):0, ("NaCl","extracto de levadura"):13, ("Fosfato","extracto de levadura"):0,
   ("Agar", "agar"):84, ("Levadura","agar"):0, ("Peptona","agar"):0, ("NaCl","agar"):5, ("Fosfato","agar"):0}

D={"fosfato": 206, "peptidos": 641, "aminoacidos":120, "extracto de levadura":145, "agar":189}

C={"fosfato":323, "peptidos": 510, "aminoacidos": 119, "extracto de levadura": 340, "agar": 826}

O=300

#-----------------
#Variables de decisión
#-----------------
for m in M:
    x={(r,m): lp.LpVariable(f"kg_del_recurso_{r}_para_el_medio_{m}",0, 300, lp.LpContinuous) for r in R for m in M}



#-----------------
#Crear el problema
#-----------------
prob=lp.LpProblem("Biocultivos", lp.LpMinimize)

#-----------------
#Restricciones
#-----------------

for r in R:
    prob+= lp.lpSum(x[r,m] for m in M)<=D[r]

for m in M:
    for n in N:
        prob+= lp.lpSum(x[r,m]*H[n,r] for r in R)>=300*J[m,n]

for m in M:
    for n in N:
        prob+=lp.lpSum(x[r,m]*H[n,r] for r in R)<=300*G[m, n]
        
for m in M:
    prob+=lp.lpSum(x[r, m] for r in R)==300


#-----------------
#Función Objetivo
#-----------------
prob+=lp.lpSum(lp.lpSum(x[r,m]*C[r]for r in R)for m in M)


#-----------------
#print(prob)
 
#Resolver el problema
prob.solve()

#Impresión de resultados
#Estado del problema
print("status: ", lp.LpStatus[prob.status])

#Valor de la función objetivo
print(" El costo es ", lp.value(prob.objective), "$")

#Valor de las variables

#-----------------
for r in R:
    for m in M: 
        print("la cantidad en kg de recurso", r, " a producir es ", lp.value(x[r, m]), " para el medio de cultivo ", m)
    


for m in M:
    cantidad = sum(lp.value(x[r,m]) for r in R)
    print("la cantidad en Kg del medio de cultivo", m, "  es ", cantidad)

for n in N:
    for m in M:
        nut= sum(H[n,r]*x[r,m])
print J    
#for m in M:
 #  for r in R:
  #      print("la cantidad del recurso ", r, " para el medio de cultivo", m, "  es ", cantidad)
