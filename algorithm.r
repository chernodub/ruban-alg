##Ядерные функции. По умолчанию используется первая.
nuclearFunction<-function(z, x = 1, r = 2, s = 3) {
  if (x == 1) {
    p<-(1-(z^r))^s
  }
  else if (x == 2) {
    p<-exp(-s*(z^r))
  }
  else if (x == 3) {
    p<-1/(s*(z^r)+1)
  }    
  else if (x == 4) {
    p<-1/(z^s)
   }
  else 
    p <- 0
  return(p)
}

f<- function(x, y)
{
  z<-function(x)
    return(-(1/((x-1)^2 + 0.2))-(1/(2*(x-2)^2 + 0.15)) -(1/(3*(x-3)^2 + 0.3)))
  return(z(x)+z(y))
}

Ruban <- function(x, deltaX, f, lower, upper, n = 500, e = 0.001, M = 1000, y = 1, q = 2, numberOfNuclearFunc = 1, r = 2, s = 100) {
  if (length(x)!= length(deltaX) || length(x)!= length(lower)  || length(x)!= length(upper))
    stop("Ошибка, не соблюдение размерности у входящих параметров (x, deltaX, upper, lower)!")
  if (n<0)
    stop("Ошибка, количество пробных точек должно быть больше 0!")
  if (e<0)
    stop("Ошибка, константа точности не должна быть меньше 0!")
  if (M<0)
    stop("Ошибка, максимальное количество итераций не должно быть меньше 0!")
  if (y<0)
    stop("Ошибка, (y) не должен быть меньше 0!")
  k<-1 ##Будем использовать для проверки первого критерия останова - превышения числа итераций
  testX<-matrix(0, n, length(x)) ## В матрицу будут записываться пробные точки. 
  ## Где n - количество пробных точек, 
  ## length(x) - количество координатых направлений
  functionValues<-rep(0, n) ## А эта матрица содержит значения функции в пробных точках
  uValues<-matrix(0, n, length(x)) ##В этой матрице будут храниться значения u для генерации пробных точек
  p<-rep(0, n) ##Ядра 
  pNorm<-rep(0, n) ##нормированные ядра
  allX<-matrix(0, M, length(x)) ##Сюда записываем все X
  allX[1,]<-x
  allResults<-f(x[1],x[2]) ##Сюда записываем все результаты. Если число циклов превысит M, то мы выберем из всех результатов наименьший
  while(k<M) {
    pastIterationX<-x ##Запоминаем текущую точку, мы будем использовать ее для условия выхода 
    plot(pastIterationX[1], pastIterationX[2], xlim=c(lower[1],upper[1]), ylim=c(lower[2],upper[2]), col="red")
    ##Шаг первый - инициализация пробных точек и расчет значений функции в них: 
    for (i in 1:n) {
      for (j in 1:length(x))
        uValues[i,j]<-runif(1, 0, 1)*2 - 1 ##Генерируем u от -1 до 1
      testX[i,]<- x + deltaX * uValues[i,]  ##Генерируем пробную точку 
      functionValues[i] <- f(testX[i,1],testX[i,2]) ##Находит значение функции в этой точке(f - какая-то функция)
      points(testX[i,1],testX[i,2], col="blue")
    }
    gmin<-rep(0, n)
    ##Шаг второй - рассчитываем ядра и нормированные ядра:
    for (i in 1:n)
      gmin[i]<-(functionValues[i] - min(functionValues)) / (max(functionValues) - min(functionValues))
    for(i in 1:n) ##Вычисляем ядра
      p[i]<-nuclearFunction(x = numberOfNuclearFunc, z = gmin[i], r=r, s=s) 
    for (i in 1:n) ##Вычисляем нормированные ядра
      pNorm[i]<-p[i]/sum(p)
    ##Шаг третий - вычисляем новый центр прямоугольника и его длину: 
    for (i in 1:length(x))
      x[i]<-x[i] + deltaX[i] * sum(sapply(1:n, function(x){uValues[x,i]*pNorm[x]}))
    print("---------------")
    cat("Результат ", k, "итерации: x=", x,"\n")
    if (x[1]<lower[1])
      x[1]<-lower[1]
    if (x[2]<lower[2])
      x[2]<-lower[2]
    if (x[1]>upper[1])
      x[1]<-upper[1]
    if (x[2]>upper[2])
      x[2]<-upper[2]
    for (i in 1:length(x))
      deltaX[i]<-y*deltaX[i] * ((sum(sapply(1:n, function(x){(abs(uValues[x,i])^(q)) * pNorm[x]})))^(1/q))
    cat("deltaX=", deltaX, "\n")
    k<-k+1 ##Увеличиваем счетчик итераций
    allX[k,]<-x
    allResults<-c(allResults,f(x[1],x[2]))
    flag <- (max(deltaX) <= e)
    if (flag == TRUE)
      break
    #if (sqrt(sum(deltaX^2)) < e) 
    #  break
  }
  cat("Сделано шагов: ", k, "\n")
  return(allX[which.min(allResults),]) 
}