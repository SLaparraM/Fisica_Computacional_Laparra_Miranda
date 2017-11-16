
__precompile__() # Este comando es para que julia precompile el paquete

module herramientas

export metodo_newton



#Calcula las raices para el pozo finito cuadrado con A=10.0

function metodo_newton(s)
    x=s;
    f(x)=x*tan(x)-sqrt(complex(10.0-x^2))+im -im;
    for i in 1:100
        df=complex(x*sec(x)^2+x/sqrt(complex(10.0-x^2))+tan(x)) #Esto se hizo para evitar errores de dominio
        x=x-f(x)/df;
    end
    return x,f(x)
end





export biseccion



#Este programa también grafica la convergencia de las raices para 4 condiciones inciales, pero utilizando el método de Bisección.

function biseccion(a,b)
    p = (a+b)/2
    f(x) = x^3-x^2-x-1
    list=zeros(50)
    for i in 1:50
        if  f(a)*f(p) > 0 
            a = p
        end;
        
        if f(b)*f(p) > 0
            b = p
        end;
        p = (a+b)/2
        list[i]=p
    end;
    return (list-list[50])./list[50]  #Este calculo se hizo para obtener el error relativo a la raiz
end






export Metodo_del_rectangulo



#Este programa calcula el valor aproximado de la integral con el Método del rectángulo, el cual consiste en sumar las áreas de todos los rectángulos formados bajo la curva. La base de cada rectángulo corresponde a una diferencia finita de un subintervalo y la altura a la función evaluada en el punto que se encuentra a la mitad del subintervalo [a,b].

function Metodo_del_rectangulo(f,c)
    i=1 
    I=0 #Se inicia la suma en I=0
    for i in 1:length(c)-1 #Se utilizó un for para ir iterando los valores de los subintervalos tomados del intervalo c
        a=c[i] # a corresponde al límite inferior del subintervalo [a,b] y c[] localiza al primer elemento del subintervalo.
        b=c[i+1] # b es el límte superior del subintervalo y c[i+1] incrementa en 1 al valor anterior y con ello localiza al                       segundo elemento del subintervalo, y así sucesivamente para todos los demás.
        I=I+(b-a)*f((a+b)/2) #La aproximación a la integral por este método, ésta dada por la fórmula I. 
    end
        return I #Devuelve el valor de la integral como la suma de todos los rectángulos de diferencia finita.
    
end





export Metodo_del_trapecio



function Metodo_del_trapecio(f,c)
    i=1
    I=0 #Se inicia la suma en I=0
    for i in 1:length(c)-1 #Se usó un for para iterar los valores del subintervalo, ya que estos valores serán usados en la suma                            de las áreas de todos los trapecios.
        a=c[i] # a es el límite superior del subintervalo 
        b=c[i+1] # b es el límite inferior del subintervalo
        I=I+(b-a)*(f(a)+f(b))/2 # I corresponde a la suma de todas las áreas de los trapecios.
    end
    return I #Devuelve el valor de la integral
end





export Metodo_Simpson



#Este progrma utliza la regla de Simpson para obtener una aproximación a la integral y sigue el mismo procedimiento que los anteriores, pero aproximando los subintervalos de f mediante polinomios de segundo grado.

function Metodo_Simpson(f,c)
    i=1
    I=0
    for i in 1:length(c)-1 #Con el uso de for se iteraron los valores del subintervalo para que la fórmula I los tome y los                                 evalue.
        a=c[i] # a corresponde al límite inferior del subintervalo
        b=c[i+1] # b es el límite superior del subintervalo
        I=I+((b-a)/6)*[f(a)+4*f((a+b)/2)+f(b)] #Consideramos el polinomio interpolador de orden dos que aproxima a la función                                                   integrando f(x) en los puntos a y b. Aplicando el método de interpolación de                                                     Lagrange, que para tres puntos interpola con un polinomio de grado 2 se obtiene                                                 la fórmula I.
    end
    return I #Devueleve el valor de I
end





export euler_explicito



function euler_explicito(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
        x = x + f(x,t)*h
        push!(listx,x) 
     end
     return listx
end





export euler_implicito



function euler_implicito(f,df,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t = i*h
       x = x -(x-f(x,t)*h)/df(x)
        push!(listx,x) 
     end
     return listx
end





export euler_implicito



function RK4(f,list,x0)
     x = x0
     h = list[2]-list[1]
     listx=[]
     push!(listx,x)
     for i in 2:length(list)
        t=i*h
        k1=f(x,t);
        k2=f(x+(h/2)*k1,t+(h/2));
        k3=f(x+(h/2)*k2,t+(h/2));
        k4=f(x+h*k3,t+h);
        x = x + (h/6)*(k1+2*k2+2*k3+k4);
        push!(listx,x) 
    end
     return listx
end

end
