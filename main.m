function varargout = main(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @main_OpeningFcn, ...
    'gui_OutputFcn',  @main_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes on button press in open_pushbutton1.
function open_pushbutton1_Callback(hObject, eventdata, handles)
global x;  %
global Fs; %
global tl;
global x2;
[filename, pathname] = uigetfile('*.wav', 'Signal');
if isequal(filename,0)
    disp('User selected Cancel')
else
    path = fullfile(pathname, filename);
    [handles.x,handles.Fs]=audioread(path);
    x=handles.x;
    Fs=handles.Fs;
    axes(handles.axes1);
    tl=[0:1/Fs:(length(handles.x)-1)/Fs];
    plot(tl,handles.x);
    title('Time Domain Waveform');
    xlabel('time/s');
    grid on;

    N=length(handles.x);
    df=Fs/N;
    w=[0:df:df*(N-1)] - Fs/2;
    X=fft(handles.x);
    X=fftshift(X);
    axes(handles.axes2);
    plot(w,abs(X)/max(abs(X)));
    axis([-10000,10000,0,1]);
    title('Frequency Domain Spectrum');
    xlabel('Frequency/Hz');
    grid on;
    x2=x;
end


% --- Executes on button press in play_pushbutton2.
function play_pushbutton2_Callback(hObject, eventdata, handles)
global x2;
global Fs;
sound(x2,Fs);
% --- Executes on button press in add_pushbutton3.
function add_pushbutton3_Callback(hObject, eventdata, handles)
global x;
global Fs;
global tl;
global x2;
axes(handles.axes1);
size(x);
t=0:1/Fs:(length(x)-1)/Fs;
Au=0.07;
fn =  get(handles.noise_edit2,'string');
fn = str2double(fn);
noise=Au*cos(2*pi*fn*t)'+0.75*Au*cos(2*pi*0.75*fn*t)'+0.5*Au*cos(2*pi*0.5*fn*t)';
x=x+noise;
plot(tl,x);
title('Time Domain Waveform');
xlabel('Time/s');
grid on;
N=length(x);
df=Fs/N;
w=[0:df:df*(N-1)] - Fs/2;
X=fft(x);
X=fftshift(X);
axes(handles.axes2);
plot(w,abs(X)/max(abs(X)));
axis([-10000,10000,0,1]);
title('Frequency Domain Spectrum');
xlabel('Frequency/Hz');
grid on;
x2=x;

% --- Executes on button press in low_pushbutton5.
function low_pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to low_pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global x;
global Fs;
global tl;
global x2;
% Copy the original audio signal data stored in variable 'x' to 'x1'.
% 'x1' will be used to store intermediate data during the filtering process to avoid directly modifying the original data 'x'.
x1=x;
% Set the passband cutoff frequency 'fp' for the low-pass filter to 3000 Hz.
fp=3000;
% Set the stopband cutoff frequency 'fs' for the low-pass filter. Here it's set to 'fp + 300', which is 3300 Hz.
fs=fp+300;
% Calculate the normalized digital frequencies 'Wp' and 'Ws' based on the sampling frequency 'Fs'.
% 'Wp' corresponds to the passband normalized frequency and 'Ws' corresponds to the stopband normalized frequency.
Wp=2*fp/Fs;
Ws=2*fs/Fs;
% Check if the normalized passband frequency 'Wp' is greater than or equal to 1.
% If so, set it to 0.99 because the range of digital frequencies should be between 0 and 1.
if(Wp>=1)
    Wp=0.99;
end
% Check if the normalized stopband frequency 'Ws' is greater than or equal to 1.
% If so, set it to 0.99 for the same reason as above.
if (Ws>=1)
    Ws=0.99;
end
% Check if the 'radiobutton1' exists in the 'handles' structure.
% This is to ensure that we can access its state property later for determining the filtering method.
if isfield(handles, 'radiobutton1')
    % Get the value of 'radiobutton1'. The value should be either 0 (unchecked) or 1 (checked).
    % This value will be used to decide which filtering method to use (Butterworth or FIR).
    radioButtonValue = get(handles.radiobutton1, 'value');
    % Check if the obtained value of 'radiobutton1' is valid (either 0 or 1).
    % If it's not valid, it indicates there might be some issues with its property settings.
    if ~(radioButtonValue == 0 || radioButtonValue == 1)
        disp('Unable to correctly obtain the state of radiobutton1. Please check its property settings.');
        % Since the state is invalid, perform a default filtering operation. Here we assume to use FIR filtering by default.
        b2=fir1(30, (Wp+Ws)/2, hamming(31));
        axes(handles.axes3);
        [h, w] = freqz(b2,1,512);
        plot(w/pi * Fs/2, 20*log(abs(h)));
        x1 = fftfilt(b2, x1);
    else
        % If the value of 'radiobutton1' is valid, perform the corresponding filtering operation based on its state.
        if radioButtonValue
            % If 'radiobutton1' is checked (value is 1), use the Butterworth filter design.
            % Calculate the order 'n' and normalized cutoff frequency 'Wn' for the Butterworth filter using the 'buttord' function.
            % The parameters 2 and 15 represent the maximum passband attenuation (2 dB) and minimum stopband attenuation (15 dB) requirements respectively.
            [n, Wn] = buttord(Wp, Ws, 2, 15);
            % Design the Butterworth low-pass filter coefficients 'b' and 'a' using the calculated parameters.
            [b, a] = butter(n, Wn, 'low');
            axes(handles.axes3);
            % Calculate the frequency response of the Butterworth filter using the 'freqz' function.
            % Get the frequency vector 'w' and the magnitude response vector 'h'.
            [h, w] = freqz(b, a);
            % Plot the magnitude response of the Butterworth filter on the specified axes.
            plot(w/pi*Fs/2, abs(h));
            % Filter the audio signal 'x1' using the designed Butterworth filter.
            x1 = filter(b, a, x1);
        else
            % If 'radiobutton1' is unchecked (value is 0), use the FIR filter design.
            % Design a FIR low-pass filter with a length of 30 and the specified cutoff frequency using the 'fir1' function.
            % The Hamming window (generated by 'hamming(31)') is used to reduce the Gibbs phenomenon.
            b2 = fir1(30, (Wp + Ws) / 2, hamming(31));
            axes(handles.axes3);
            [h, w] = freqz(b2, 1, 512);
            plot(w/pi*Fs/2, 20*log(abs(h)));
            x1 = fftfilt(b2, x1);
        end
    end
else
    disp('radiobutton1 does not exist. Please check the GUI design.');
    % As 'radiobutton1' doesn't exist, perform a default FIR filtering operation.
    b2=fir1(30, (Wp+Ws) / 2, hamming(31));
    axes(handles.axes3);
    [h, w] = freqz(b2, 1, 512);
    plot(w/pi*Fs/2, 20*log(abs(h)));
    x1 = fftfilt(b2, x1);
end
% Plot the time-domain waveform of the filtered audio signal 'x1' on the specified axes.
% Set the title and x-axis label for the time-domain waveform plot.
axes(handles.axes5);
plot(tl, x1);
title('Time Domain Waveform');
xlabel('Time/s');
% Calculate the frequency axis data for plotting the frequency-domain spectrum.
% Get the length of the filtered audio signal 'x1' and calculate the frequency resolution 'df'.
% Then create the frequency axis vector 'w'.
N=length(x1);
df=Fs/N;
w = [0:df:df*(N-1)]-Fs/2;
% Perform the Fast Fourier Transform (FFT) on the filtered audio signal 'x1'.
X = fft(x1);
% Shift the zero-frequency component to the center of the spectrum for better visualization.
X = fftshift(X);
% Plot the frequency-domain spectrum of the filtered audio signal 'x1' on the specified axes.
% Normalize the magnitude of the spectrum by dividing it by the maximum magnitude.
% Set the axis limits, title and x-axis label for the frequency-domain spectrum plot and turn on the grid.
axes(handles.axes6);
plot(w, abs(X) / max(abs(X)));
axis([-10000, 10000, 0, 1]);
title('Frequency Domain Spectrum');
xlabel('Frequency/Hz');
grid on;
x2 = x1;


% --- Executes during object creation, after setting all properties.
function low_pushbutton5_CreateFcn(hObject, eventdata, handles)
% --- Executes on button press in pushbutton6.


function noise_edit2_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of noise_edit2 as text


% --- Executes during object creation, after setting all properties.
function noise_edit2_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function add_pushbutton3_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in guass_pushbutton8.
function guass_pushbutton8_Callback(hObject, eventdata, handles)

global x;
global Fs;
global tl;
global x2;
axes(handles.axes1);
x = awgn(x,15);
plot(tl,x);
title('Time Domain Waveform');
xlabel('Time/s');
grid on;
N=length(x);
df=Fs/N;
w=[0:df:df*(N-1)] - Fs/2;
X=fft(x);
X=fftshift(X);
axes(handles.axes2);
plot(w,abs(X)/max(abs(X)));
axis([-10000,10000,0,1]);
title('Frequency Domain Spectrum');
xlabel('Frequency/Hz');
grid on;
x2=x;


function high_pushbutton9_Callback(hObject, eventdata, handles)

global x;
global Fs;
global tl;
global x2;

x1=x;
fp = 3000;
fs = fp-300;
Wp = 2*fp/Fs;
Ws = 2*fs/Fs;
if(Wp >= 1)
    Wp = 0.99;
end
if(Ws >= 1)
    Ws = 0.99;
end
if get(handles.radiobutton1,'value')
    [n, Wn]=buttord(Wp,Ws, 2, 15);
    [b, a]=butter(n, Wn,'high');
    axes(handles.axes3);
    [h,w]=freqz(b,a);
    plot(w/pi*Fs/2,abs(h));
    x1=filter(b,a,x1);
else
    b2=fir1(30, (Wp+Ws)/2,'high',hamming(31));
    axes(handles.axes3);
    [h,w]=freqz(b2, 1,512);
    plot(w/pi*Fs/2,20*log(abs(h)));
    x1=fftfilt(b2,x1);
end;
axes(handles.axes5);
plot(tl,x1);
title('Time Domain Waveform');
xlabel('Time/s');
N=length(x1);
df=Fs/N;
w=[0:df:df*(N-1)] - Fs/2;
X=fft(x1);
X=fftshift(X);
axes(handles.axes6);
plot(w,abs(X)/max(abs(X)));
axis([-10000,10000,0,1]);
title('Frequency Domain Spectrum');
xlabel('Frequency/Hz');
grid on;
x2=x1;

% --- Executes on button press in bandpass_pushbutton10.
function bandpass_pushbutton10_Callback(~, eventdata, handles)

global x;
global Fs;
global tl;
global x2;

x1=x;
fp = [1500,3000];
fs = [1000,3300];
Wp = 2*fp/Fs;
Ws = 2*fs/Fs;
if(Wp >= 1)
    Wp = 0.99;
end
if(Ws >= 1)
    Ws = 0.99;
end
if get(handles.radiobutton1,'value')
    [n, Wn]=buttord(Wp,Ws, 2, 15);
    [b, a]=butter(n, Wn,'bandpass');
    axes(handles.axes3);
    [h,w]=freqz(b,a);
    plot(w/pi*Fs/2,abs(h));
    x1=filter(b,a,x1);
else
    b2=fir1(30, (Wp+Ws)/2,hamming(31));
    axes(handles.axes3);
    [h,w]=freqz(b2, 1,512);
    plot(w/pi*Fs/2,20*log(abs(h)));
    x1=fftfilt(b2,x1);
end;
axes(handles.axes5);
plot(tl,x1);
title('Time Domain Waveform');
xlabel('Time/s');
N=length(x1);
df=Fs/N;
w=[0:df:df*(N-1)] - Fs/2;
X=fft(x1);
X=fftshift(X);
axes(handles.axes6);
plot(w,abs(X)/max(abs(X)));
axis([-10000,10000,0,1]);
title('Frequency Domain Spectrum');
xlabel('Frequency/Hz');
grid on;
x2=x1;

% --- Executes on button press in stop_pushbutton11.
function stop_pushbutton11_Callback(hObject, eventdata, handles)
global x;
global Fs;
global tl;
global x2;

x1=x;
fp = [1200,3300];
fs = [1500,3000];
Wp = 2*fp/Fs;
Ws = 2*fs/Fs;
if(Wp >= 1)
    Wp = 0.99;
end
if(Ws >= 1)
    Ws = 0.99;
end
if get(handles.radiobutton1,'value')
    [n, Wn]=buttord(Wp,Ws, 2, 15);
    [b, a]=butter(n, Wn,'stop');
    axes(handles.axes3);
    [h,w]=freqz(b,a);
    plot(w/pi*Fs/2,abs(h));
    x1=filter(b,a,x1);
else
    b2=fir1(30, (Wp+Ws)/2,'stop',hamming(31));
    axes(handles.axes3);
    [h,w]=freqz(b2, 1,512);
    plot(w/pi*Fs/2,20*log(abs(h)));
    x1=fftfilt(b2,x1);
end;
axes(handles.axes5);
plot(tl,x1);
title('Time Domain Waveform');
xlabel('Time/s');
N=length(x1);
df=Fs/N;
w=[0:df:df*(N-1)] - Fs/2;
X=fft(x1);
X=fftshift(X);
axes(handles.axes6);
plot(w,abs(X)/max(abs(X)));
axis([-10000,10000,0,1]);
title('Frequency Domain Spectrum');
xlabel('Frequency/Hz');
grid on;
x2=x1;


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in invert_pushbutton12.
function invert_pushbutton12_Callback(hObject, eventdata, handles)
global x;
global x2;

N=length(x);
for i=1:1:N
    x2(i)= x(N-i+1);
end
sound(x2);


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
close;

% --- Executes on button press in radiobutton2.
function radiobutton5_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton1.
function radiobutton4_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of radiobutton1
