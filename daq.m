s = daq.createSession('ni');
addAnalogInputChannel(s,'P6363','ai8','Voltage');
lh = addlistener(s,'DataAvailable',@plotData); 
s.IsContinuous = true;

startBackground(s);

%%
stop(s)
%%
function plotData(src,event)
     plot(event.TimeStamps,event.Data)
end