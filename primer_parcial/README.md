# Práctica 1: Valores y vectores propios

El repositorio contiene los seis ejercicios indicados en la práctica 1 de
Cálculo Numérico en el ITAM (otoño 2017).

Para el proyecto utilizamos Julia, con sólo dos librerías adicionales a la base:
`Distritbutions` en el ejercicio 6 y `MatrixDepot` en el ejercicio 5. En ambos
casos, la primera línea de cada ejercicio (que dejamos comentada) instala la
librería adecuada.

Fuera del contenido estándar de la clase, queremos recalcar las siguientes
particularidades de nuestro proyecto:

-El método de la potencia-Rayleigh tiene como condición de finalización
alternativa a superar el número de iteraciones el que la matriz $A-\rho I$,
donde $\rho$ es el shift de Rayleigh actual, sea singular.

-Para los métodos QR, usamos factorización de Hessenberg.

-Para el método QR dinámico, usamos shift de Wilkinson.

Todas las funciones auxiliares utilizadas están en la carpeta lib, y están
documentadas en Julia. Para acceder a su documentación, puede usarse
`? func` desde la terminal de Julia, donde `func` es el nombre de la función.
Cada ejercicio tiene comentarios destinados a que los ejercicios sean
autocontenidos, para evitar en la medida de lo posible cambiar entre el
código y el trabajo escrito.
