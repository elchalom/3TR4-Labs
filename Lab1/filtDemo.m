%
% Filter Demo
%
%

clear all
clear functions
close all

load FiltSig


%default parms
Mtitle = 'Filter Demonstration';
Stitle = 'Original Signal Spectrum';
Ftitle = 'Filter Response';
Rtitle = 'Resulting Response';

Mposition=[90 520 345 180];
Sposition=[30 60  430  315];
Fposition=[500 440 430 315];
Rposition=[500 60  430 315];

%open windows

SfigNumber=figure( ...
        'Name',Stitle, ...
        'NumberTitle','off', ...
        'NextPlot','add', ...
	'Visible','off', ...
	'Position',Sposition, ...
	'Colormap',[]);
plot( spFreq, Pspeech  );xlabel('Frequency');
ylabel('Power Spectrum Magnitude (dB)');grid
set(SfigNumber,'NextPlot','new');


playSpeech = 'soundsc( orgspeech)';


FfigNumber=figure( ...
        'Name',Ftitle, ...
        'NumberTitle','off', ...
        'NextPlot','new', ...
	'Visible','off', ...
	'Position',Fposition, ...
	'Colormap',[]);


RfigNumber=figure( ...
        'Name',Rtitle, ...
        'NumberTitle','off', ...
        'NextPlot','new', ...
	'Visible','off', ...
	'Position',Rposition, ...
	'Colormap',[]);


MfigNumber=figure( ...
        'Name',Mtitle, ...
        'NumberTitle','off', ...
        'NextPlot','new', ...
	'Visible','off', ...
	'Position',Mposition, ...
	'Colormap',[]);


%  Configure layout of Window
 top=0.95;   left=0.05;  right=0.95;  bottom=0.05;
 labelHt=0.05;    spacing=0.005;

% Create the Text Window frame
    frmBorder=0.02;
    frmPos=[left-frmBorder bottom-frmBorder ...
        (right-left)+2*frmBorder (top-bottom)+2*frmBorder];



%====================================
% Information for all buttons

labelColor=[0.8 0.8 0.8];
top=0.95; bottom=0.1; yInitPos=0.05;  btnWid=0.4;  btnHt=0.20;

% Spacing between the button and the next command's label
spacing=0.035;


buttonPos = [ yInitPos top-btnHt btnWid btnHt];
%====================================
%Button commands
nullCmd = ' ';
lPassCmd =  ['exButton(1,fbutton);'...
		'figure(SfigNumber); set(gcf,''NextPlot'',''add'');'...
		'plot( spFreq, Pspeech  ); xlabel(''Frequency'');' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');grid; ' ...
		'set(SfigNumber,''NextPlot'',''new'');figure(FfigNumber);' ...
                'set(gcf,''NextPlot'',''add''); ' ...
		'plot(spFreq,lfP); xlabel(''Frequency'');yMin(-70);' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');' ...
		'grid; set(gcf,''NextPlot'',''new'');drawnow;' ...
		'figure(RfigNumber);' ...
                'set(gcf,''NextPlot'',''add'') ;' ...
		'plot(spFreq,lpP);xlabel(''Frequency'');yMin(-90);' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');' ...
		'grid;drawnow;set(gcf,''NextPlot'',''new'');' ...
        'soundsc(speech); pause(3);' ...
		'soundsc(lpspeech)' ...
            ];
bPassCmd =  ['exButton(2,fbutton);'...
		'figure(SfigNumber); set(gcf,''NextPlot'',''add'');'...
		'plot( spFreq, Pspeech  ); xlabel(''Frequency'');' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');grid; ' ...
		'set(SfigNumber,''NextPlot'',''new'');figure(FfigNumber);' ...
                'set(gcf,''NextPlot'',''add''); ' ...
		'plot(spFreq,bfP); xlabel(''Frequency'');yMin(-70);' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');' ...
		'grid; drawnow;set(gcf,''NextPlot'',''new'');' ...
		'figure(RfigNumber);' ...
                'set(gcf,''NextPlot'',''add'') ;' ...
		'plot(spFreq,bpP);xlabel(''Frequency'');yMin(-90);' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');' ...
		'grid;drawnow;set(gcf,''NextPlot'',''new'');' ...
        'soundsc(speech); pause(3);' ...
		'soundsc( bpspeech)' ...
            ];
hPassCmd =  ['exButton(3,fbutton);'...
		'figure(SfigNumber); set(gcf,''NextPlot'',''add'');'...
		'plot( spFreq, Pspeech  ); xlabel(''Frequency'');' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');grid; ' ...
		'set(SfigNumber,''NextPlot'',''new'');' ...
		'figure(FfigNumber);' ...
                'set(gcf,''NextPlot'',''add''); ' ...
		'plot(spFreq,hfP); xlabel(''Frequency'');yMin(-70);' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');' ...
		'grid; drawnow;set(gcf,''NextPlot'',''new'');' ...
		'figure(RfigNumber); set(gcf,''NextPlot'',''add'') ;' ...
		'plot(spFreq,hpP);xlabel(''Frequency'');yMin(-90);' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');' ...
       'grid;drawnow;set(gcf,''NextPlot'',''new'');' ...
       'soundsc(speech); pause(3);' ...
		'soundsc(hpspeech)' ...
            ];
nPassCmd =  ['exButton(4,fbutton);'...
		'figure(SfigNumber); set(gcf,''NextPlot'',''add'');'...
		'plot( spFreq, ntPspeech  ); xlabel(''Frequency'');' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');grid; ' ...
		'set(SfigNumber,''NextPlot'',''new'');' ...
 		'figure(FfigNumber);  set(gcf,''NextPlot'',''add''); ' ...
		'plot(spFreq,nfP); xlabel(''Frequency'');yMin(-70);' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');' ...
		'grid; drawnow;set(gcf,''NextPlot'',''new'');' ...
		'figure(RfigNumber);' ...
                'set(gcf,''NextPlot'',''add'') ;' ...
		'plot(spFreq,npP);xlabel(''Frequency'');yMin(-90);' ...
		'ylabel(''Power Spectrum Magnitude (dB)'');' ...
		'grid;drawnow;set(gcf,''NextPlot'',''new'');' ...
	   'soundsc(ntspeech); pause(3);' ...
       'soundsc(npspeech)' ...
            ];

genIRQ = ' close all';


%====================================
% The CONSOLE frame
frmBorder=0.02;
uicontrol( ...
        'Style','frame', ...
        'Units','normalized', ...
        'Position',frmPos, ...
	'BackgroundColor',[0.5 0.5 0.5]);

fbutton(1) = uicontrol( ...
	'Style', 'pushbutton', ...
        'Units','normalized', ...
        'Position',buttonPos, ...
	'String','Low Pass', ...
	'Value',0, ...
	'CallBack', lPassCmd, ...
	'BackgroundColor',[0.7 0.7 0.7]);

buttonPos(2)=buttonPos(2)-(btnHt+spacing);


fbutton(2) = uicontrol( ...
	'Style', 'pushbutton', ...
        'Units','normalized', ...
        'Position',buttonPos, ...
	'String','Band Pass', ...
	'Value',0, ...
	'CallBack', bPassCmd, ...
	'BackgroundColor',[0.7 0.7 0.7]);

buttonPos(2)=buttonPos(2)-(btnHt+spacing);


fbutton(3) = uicontrol( ...
	'Style', 'pushbutton', ...
        'Units','normalized', ...
        'Position',buttonPos, ...
	'String','High Pass', ...
	'Value',0, ...
	'CallBack', hPassCmd, ...
	'BackgroundColor',[0.7 0.7 0.7]);

buttonPos(2)=buttonPos(2)-(btnHt+spacing);


fbutton(4) = uicontrol( ...
	'Style', 'pushbutton', ...
        'Units','normalized', ...
        'Position',buttonPos, ...
	'String','Notch', ...
	'Value',0, ...
	'CallBack', nPassCmd, ...
	'BackgroundColor',[0.7 0.7 0.7]);


 
buttonPos(1)=buttonPos(1)+(btnWid+2*spacing);

    uicontrol( ...
        'Style','pushbutton', ...
        'Units','normalized', ...
        'Position',buttonPos, ...
        'String','exit', ...
        'Callback',genIRQ);


    uicontrol( ...
        'Style','pushbutton', ...
        'Units','normalized', ...
        'Position',[yInitPos+(btnWid+2*spacing) top-btnHt btnWid btnHt], ...
        'String','Original Speech', ...
        'Callback','soundsc( speech)');


%display default curves
%eval(lPassCmd);

%activate all display figures
figure(MfigNumber);
set(MfigNumber,'Visible','on');
set(SfigNumber,'Visible','on');
set(FfigNumber,'Visible','on');
set(RfigNumber,'Visible','on');

    
