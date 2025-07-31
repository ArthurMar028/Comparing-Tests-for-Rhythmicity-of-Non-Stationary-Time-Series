#Utilizando um dado real "Posicao_vs_tempo.csv", responda as seguintes perguntas, criando um algoritmo otimizado para a análise deste tipo de dado por meio de diferentes métodos:

#Qual o período dominante das oscilações em taxa de crescimento nestes dados de posição (em pixels) vs tempo (em minutos)? 
#Qual a taxa média e máxima de crescimento em micrômetros por segundo, sendo que 1 pixel é 0.1084 micrômetros?

#Referências: Damineli et al. (2017), Leise (2013), Leise (2017), Rösch & Schmidbauer (2016).

###########################################################################

#0. Data loader

data <- read.csv("Posicao_vs_tempo.csv") #As it is a position x time series, there is already a oscilatory component, but with a linear (I think) trend
head(data)
posicao <- data$Posicao
time <- data$Tempo #minutos
plot(time, posicao, type="l", col="red3")

#A) Preprocessing Data

#1. Trend removal

#1.1. Making a low-degree polynomial fit to the ts

#linear

linear_model <- loess(posicao ~ time, degree = 1, span = 0.1)
plot(time, linear_model$residuals, type="l", col="red3") #detrended series

#linear_model <- lm(y ~ time)
#y_linear_trend <- predict(linear_model) #I need to find a way of making this a list without the numbers as labels (maybe it is related to lm() or to the data analysed itself)
#abline(linear_model, col="blue")


#stat_l <- summary(linear_model) #To compare R values?
#summary(linear_model)

#quadratic

quad_model <- loess(posicao ~ time, degree = 2, span = 0.1)
plot(time, quad_model$residuals, type="l", col="blue") #detrended series by a 2nd degree polynomial


#quad_model <- lm(y ~ time + I(time^2))
#curve(predict(quadratic_model, newdata = data.frame(x = x)), add = TRUE, col = "green4") #plot can deal with only 2 parameters (reason for using curve); predict is used to transform a matrix (lm)into a function 

#summary(quadratic_model) #both R values are high, but I think the linear fit makes more sense as the values that makes it more "imprecise" are those in the beggining and end of the series (? or i might be wrong)

#y_quad_model <- predict(quad_model)

##(not necessary when using residuals) 1.2. Subtracting the treand (might be a running average (idk how to do so) or a low-degree polynomial fit)

#detrended_linear <- y - y_linear_trend #? weird result - I think only a wavelet analysis would make it fit more adequately to be around 0

#detrended_quad <- y - y_quad_model

#1.2. plotting the growth rate (the difference of the position of the residuals over time - which mean, the difference of a fit and the actual curve)

par(mfrow=c(1,3))

plot(time[-1], (diff(posicao)/diff(time)[1])*0.1084, type="l", col="blue3") #Calcular todas as vezes diff(time)[1] faz com que ele fique recalculando (CUSTOSO E NÃO É boa prática)
]]
abline(h=median((diff(posicao)/diff(time))*0.1084))
plot(time[-1], (diff(linear_model$residuals)/diff(time)[1])*0.1084, type="l", col="red3")
abline(h=median((diff(linear_model$residuals)/diff(time)[1])*0.1084))
plot(time[-1], (diff(quad_model$residuals)/diff(time)[1])*0.1084, type="l", col="green3")
abline(h=median((diff(quad_model$residuals)/diff(time)[1])*0.1084))

par(mfrow=c(1,1))

#2. Outlier removal (?) - beggining and end? - necessary?

#lin_out <- linear_model$residuals[152,380]

#Plot: "Dashboard" with (original series + series with peak and ISI + Fourier with main peak + ACF with first peak after lag 0)

#Adding a tendency

#3) 
#a) Testar outras funções do R para fazer a estimação da periodicidade como 'spectrum' (comparar spec.pgram & spec.ar) e entender como extrair somente os picos significativos. 
#b) Testar implementações de Short Time Fourier Transforms (como no pacote signal e função specgram) para o caso de séries não estacionárias.
#c) Utilizar o pacote WaveletComp para obter estimativas temporalmente explícitas de séries não estacionárias. Esse guia é excelente, recomendo seguir e utilizar ainda mais pra frente nas análises:
  #http://www.hs-stat.com/projects/WaveletComp/WaveletComp_guided_tour.pdf

Fourier_transform <- function(time = time, signal = diff(posicao), colour = "red3"){
  
  #Fourier transform

  fourier_module <- Mod(
    fft(signal)) #As I'm looking for the growth rate oscillations, signal should already be diff(position)
  
  time_step <- diff(time)[1]
  Fs <- 1/time_step
  n <- length(time)
  freq_bins <- (0:
                  (n-1)) * Fs/n
  
  if(n %% 2 == 0) {
    metade <- n/2 + 1
  } else {
    metade <- (n + 1)/2
  }
  
  ft <- cbind.data.frame(freq_bins = freq_bins[1:metade],
                         amp = fourier_module[1:metade])
  
  ft$amp[1] <- NA #0 Hz
  
  #Finding per_fourier
  
  per <- 1/freq_bins[1:metade] 
  main_per_fourier <- per[which.max(ft$amp)]
  
  plot(per, ft$amp, 
       type = "h", col = colour,
       xlim = c(0, 2*main_per_fourier),
       sub = paste("O período principal segundo a Transformada de Fourier é",
                                                      round(main_per_fourier, digits = 3), "min")) #How could I take the mos significant ones (not only the main?)
  
}

Autocorrel <- function(time = time, signal = diff(posicao), colour = "red3"){
  
  autocorrel <- acf(signal, lag.max = length(signal), plot=FALSE)
  peaks <- diff(sign(diff(autocorrel$acf)))
  position_peaks <- (which(peaks==-2)+1) #desconsidera o primeiro pico (lag=0) e é só "+1" porque trabalho com a diff do signal (growth rate)
  time_peaks <- time[position_peaks]
  max_1st <- time_peaks[1]
  avg_acf <- mean(diff(time_peaks))
  plot(autocorrel, 
       col = colour,
       sub = paste("O período principal segundo o primeiro pico de autocorrelação após lag[1] é",
                   round(max_1st, digits = 3), "min"))
  
}

ISI <- function(time = time, signal = diff(posicao), colour = "red3"){
  
  dif_2nd <- diff(sign(diff(signal)))
  peaks_position <- which(dif_2nd==-2) + 1
  time_peaks <- time[peaks_position]
  avg_isi <- mean(diff(time_peaks))
  
  if(length(time)!=length(signal)){
    
    plot(time[-1], signal,
         type="l", col=colour,
         sub = paste("O período estimado a partir da 'Inter-Spike Interval' média é", round(avg_isi, digits=3), "min"))
    
  } else { 
  
  plot(time, signal,
       type="l", col = colour,
    sub = paste("O período estimado a partir da 'Inter-Spike Interval' média é", round(avg_isi, digits=3), "min"))
  }
  abline(v=time_peaks[1:(length(time_peaks)/5)], col="purple", lwd=0.01) #It's not exactly in the peaks
  
}

#Specgram? "Aí no spec.pgram tem uma opção de suavizar o periodograma e aí tem uma porção de métodos possíveis"

#Wavelet comp



  
#4) FINALMENTE analisar os dados utilizando as funções desenvolvidas até aqui e aprimorando elas conforme o necessário. Se precisar analisar sincronização o guia do WaveletComp é perfeito pra aprender.

Plot_growthvstime <- function(time = time, signal = posicao, colour="red3", median = TRUE, fitted = FALSE){
  
  growth_rate <- diff(signal)/mean(diff(time))*0.1084 #(pixel/min)*(micron/pixel)
  
  if(median==TRUE & fitted==FALSE){
    plot(time[-1], growth_rate, 
       type="l", col=colour,
       ylab = "Taxa de Crescimento (micron*min^-1)", xlab = "Tempo (min)", 
       main = paste("A mediana da taxa de crescimento é de", round(median(growth_rate), digits = 3), "microns*min^-1"), 
       sub = paste("A taxa de crescimento mínima é de", round(min(growth_rate), digits = 3), "e a máxima", round(max(growth_rate), digits = 3), "microns*min^-1"))
    abline(h=median(growth_rate), col=colour)
  } else if (median==FALSE & fitted==FALSE){
  plot(time[-1], growth_rate, 
       type="l", col=colour,
       ylab = "Taxa de Crescimento (micron*min^-1)", xlab = "Tempo (min)")
  } else {
    plot(time[-1], growth_rate, 
         type="l", col=colour,
         ylab = "Desvio da taxa de Crescimento média (micron*min^-1)", xlab = "Tempo (min)")
  }
  
}

#Final plot

#taxa média e máxima de crescimento em micrômetros por segundo, sendo que 1 pixel é 0.1084 micrômetros (abaixo dos gráficos de growth rate vs. time)

Plot_growthvstime(time, posicao) #Growth vs. time

par(mfrow = c(4, 3), oma = c(0, 0, 0, 0))

#fazer em funções para ficar retestando mais facilmente! - comparar raw, linear, quad e 3rd? (ou a função de posição, não sei)

# Plot Taxas de crescimento e desvios
Plot_growthvstime(time, posicao
                  )
Plot_growthvstime(time, linear_model$residuals, 
                  colour = "blue3", median = FALSE, fitted = TRUE)
mtext("Retirada tendência a partir de fit linear", 
      side = 3, line = 1, adj = 0.5, cex = 1.2, font = 1)
Plot_growthvstime(time, quad_model$residuals, 
                  colour = "green4", median = FALSE, fitted = TRUE)
mtext("Retirada tendência a partir de fit quadrática", 
      side = 3, line = 1, adj = 0.5, cex = 1.2, font = 1)

# Plot Transformadas de Fourier
Fourier_transform(time, diff(posicao)
                  )
Fourier_transform(time, diff(linear_model$residuals), 
                  colour = "blue3")
Fourier_transform(time, diff(quad_model$residuals), 
                  colour = "green4")

# Plot acf
Autocorrel(time, diff(posicao)
           )
Autocorrel(time, diff(linear_model$residuals), 
           colour = "blue3")
Autocorrel(time, diff(quad_model$residuals), 
           colour = "green4")


# Plot ISI
ISI(time, diff(posicao)
    )
ISI(time, diff(linear_model$residuals), 
    colour = "blue3")
ISI(time, diff(quad_model$residuals), 
    colour = "green4")

par(mfrow = c(1, 1))




