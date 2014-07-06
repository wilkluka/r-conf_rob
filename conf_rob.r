# Author: Łukasz Wilk
# Available as Freeware
# All rights reserved
# if you wish to use this code contact first:  
#    lukasz.s.wilk a#t gmail.com 

gir <- function( 
	n = 25, 
	jumpP = 0.001, 
	alpha = 0.05, 
	intervals = c("cp", "wald", "wilson", "midpv"), 
	style = c(1, 2, 3, 4), 
	color = c(1, 2, 3, 4) 
	) {
	
	#function for plotnig computed solution
	rys <- function( input, style, color ) {
		
		intervals = input$intervals
		alpha = input$alpha
		x = input$x
		
		plot( c(0, 1), c(0.7, 1), type = "n", xlab = "probability of success \"p\"",
		ylab = "probability that interval includes p" )
		if ( sum( "cp" == intervals ) == 1 ) { #check if user want to draw this plot
			p = sum( ("cp" == intervals) * c( 1:length(intervals) ) ) #this is for having same order in result as in intervals
			lines( x, input$inputCP, type = "l", lty = style[p], col = color[p] )
		}
		if( sum( "wald" == intervals ) == 1 ) {
			p = sum( ("wald" == intervals) * c( 1:length(intervals) ) )
			lines( x, input$inputWald, type = "l", lty = style[p], col = color[p] )
		}
		if( sum( "wilson" == intervals ) == 1 ) {
			p=sum( ("wilson" == intervals) * c( 1:length(intervals) ) )
			lines( x, input$inputWilson, type = "l", lty = style[p], col = color[p] )
		}
		if( sum( "midpv" == intervals ) == 1 ) {
			p=sum( ("midpv" == intervals) * c( 1:length(intervals) ) )
			lines( x, input$inputMidPV, type = "l", lty = style[p], col = color[p])
		}
		intervals[intervals == "cp"]="Clopper-Pearson"
		intervals[intervals == "wald"]="Wald"
		intervals[intervals == "wilson"]="Wilson"
		intervals[intervals == "midpv"]="mid-P-value"
		#here we draw line 1-alpha
		lines( c(0, 1), c(1-alpha, 1-alpha), type = "l", lwd = 1, col = 1 )
		legend(0.4, 0.8, intervals, cex = 0.8, lty = style, title = "types of confidence intervals", col = color)
		title(paste("Comparison of confidence intervals robustness for n =",n))
	}


  #function for computing roustnes of intervals
	gen <- function( n, jumpP, alpha, intervals ) {

		left.end <- function(p, x, n, alpha) {
			pbeta(p, x, n-x+1) - alpha/2 - dbinom(x, n, p)/2
		}
		
		right.end <- function(p, x, n, alpha){
			1-pbeta(p, x+1, n-x) - alpha/2 - dbinom(x, n, p)/2
		}
		
		library(Hmisc)
		q = 1/jumpP
		inputCP = c(0:q)
		inputWald = c(0:q)
		inputWilson = c(0:q)
		inputMidPV = c(0:q)
		#  these are vectors with "accuracy"
		
		# here we are starting to count intervals 
		intervalCP = binconf( c(0:n), n, method = "exact", alpha = alpha)
		intervalWald = binconf( c(0:n), n, method="asymptotic", alpha = alpha)
		intervalWilson = binconf( c(0:n), n, method="wilson", alpha = alpha)
		intervalMPV = 0
		for(k in 1:n) {
			intervalMPV = c(intervalMPV, uniroot(left.end, interval = c(0,k/n), x = k, n = n, alpha = alpha)$root )
		}
		
		for( k in 0:(n-1) ) {
			intervalMPV = c(intervalMPV, uniroot(right.end, interval = c(k/n,1), x = k, n = n, alpha = alpha)$root )
		}
		intervalMPV = c(intervalMPV, 1)
		
		if( sum( "cp" == intervals ) > 0 ) {
			p = 0
			for( i in 0:(q/2) ) {
				p = jumpP * i
				# p - current "x-axis"
				# q - amount of p-steps
				densityB = dbinom(c(0:n), n, p)
				inputCP[q+1-i] = inputCP[i+1] = sum( (intervalCP[ c( (n+2):(2*n+2) ) ]<=p) * (intervalCP[ c( (2*n+3):(3*n+3) ) ]>=p) * densityB )
			}
		}
		if( sum( "wald" == intervals ) > 0 ) {
			p = 0
			for( i in 0:(q/2) ) {
				p = jumpP * i
				densityB = dbinom(c(0:n), n, p)
				inputWald[q+1-i] = inputWald[i+1] = sum( (intervalWald[ c( (n+2):(2*n+2) ) ] < p) 
				        * (intervalWald[ c( (2*n+3):(3*n+3) ) ]>p) * densityB )
			}
		}
		if( sum( "wilson" == intervals ) > 0 ) {
			p = 0
			for( i in 0:(q/2) ) {
				p = jumpP * i
				densityB = dbinom(c(0:n), n, p)
				inputWilson[q+1-i] = inputWilson[i+1] = sum( (intervalWilson[ c( (n+2):(2*n+2) ) ]<p) * (intervalWilson[ c( (2*n+3):(3*n+3) ) ]>p) * densityB )
			}
		}
		if( sum( "midpv" == intervals ) > 0 ) {
			p = 0
			for( i in 0:(q/2) ) {
				p = jumpP * i
				densityB = dbinom(c(0:n), n, p)
				inputMidPV[q+1-i] = inputMidPV[i+1] = sum( (intervalMPV[ c(1:(n+1) ) ]<p) * (intervalMPV[ c( (n+2):(2*n+2) ) ]>p) * densityB )
			}
		}
			# assembling results
		result <- list(
			n = n, 
			jumpP = jumpP, 
			x = seq(0, 1, length = q+1), 
			inputCP = inputCP, 
			inputWilson = inputWilson, 
			inputWald = inputWald, 
			inputMidPV = inputMidPV,
			alpha = alpha,
			intervals = intervals
		)

		return(result)
	}

	rys( gen(n, jumpP, alpha, intervals), style, color)

#input (in function rys) - heve we input gen function output
#n - sample vector size
#jumpP - plot (and solution) resolution
#alpha - same as alpha elswhere (+ it draws line at 1-alpha)
#intervals - vector which methods to use "cp" "wald" "wilson" "midpv"
#style vektor lty(for plot)
#color vektor col(for plot)
}
gir()

# All rights reserved
