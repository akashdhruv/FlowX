from flowx import __environment__

if __environment__ is 'serial': from flowx.domain.serial import domain_main
if __environment__ is 'parallel': from flowx.domain.parallel import domain_main  
