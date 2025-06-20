function progress_msg(varargin)
%PROGRESS_MSG Progress messages and timing.
% 
% PROGRESS_MSG('msg',I,N,opt) where I = 0 or -1 initilises the messages. 
%   'msg'   [optional] string to print
%   I       [optional] if I=0 initialises and prints '0%', '0/N', or '0'
%                      if I=-1 suppresses the initial progress message and
%                      ignores options of the immediately following 
%                      initialisation message                      
%   N       [optional] if N=0 prints 'I', if N>0 prints 'I/N' 
%                      if absent, prints 'I%' (default)
%   opt     [optional] structure with the following fields and defaults:
%      .log_level = 2        controls verbosity, see EIDORD_MSG for details
%      .final_msg            customises the text of the final message
%                   'Done'   default
%                   ''       last message will not be erased, but time will
%                            be printed
%                   'none'   completely suppress the final message
%      .num_digits = 6       controls the text field width for N=0
%
% PROGRESS_MSG(I) displays I%
% PROGRESS_MSG(I,N) displays I/N
% PROGRESS_MSG(I,0) displays I
%
% PROGRESS_MSG(Inf) prints the final message, including the time elapsed
%  since initialisation, measured using tic/toc.
% PROGRESS_MSG('msg',Inf) prints the 'msg' instead of the opt.final_msg
%
% Example:
%   progress_msg('The big loop:', 0,10);         % prints 0/10
%   for i = 1:10
%       progress_msg(i,10);     % e.g.   3/10
%       pause(.2); % think
%   end
%   progress_msg(Inf)           % final message, e.g. Done (2.015 s)
%
%   opt.final_msg = 'none'                   % suppress final message 
%   progress_msg('Start calculations...',-1);% suppress next initialisation
%   res = function_using_progress_msg;
%   progress_msg(sprintf('Size is %d',numel(res)), Inf);
%
% SEE ALSO EIDORS_MSG, TIC, TOC

% (C) 2015 Bartlomiej Grychtol - all rights reserved Swisstom AG
% License: GPL version 2 or 3
% $Id: progress_msg.m 6981 2024-11-12 17:34:27Z aadler $

% >> SWISSTOM CONTRIBUTION <<

t0 = tic;

if nargin > 0 && ischar(varargin{1}) && strcmp(varargin{1},'UNIT_TEST')
    varargin(1) = [];
    do_unit_test(varargin{:});
    return
end

persistent pvar;

nargs = nargin;
opt = default_opts;
if nargs > 0 && isstruct(varargin{end})
    tmp = varargin{end};
    try
        if tmp.log_level > eidors_msg('log_level')
            return % don't even process options 
        end
    end
    fld = fieldnames(tmp);
    for i = 1:numel(fld)
        opt.(fld{i}) = tmp.(fld{i});
    end
    nargs = nargs-1;
    varargin(end) = [];
end


msg = '';
if nargs > 0 && ischar(varargin{1})
    if nargs > 1 && isinf(varargin{2})
        msg = varargin{1}; % custom final message from function call
    else
        msg = [repmat(' ',1,opt.log_level) datestr(now) ': ' varargin{1}]; % Fixed text
    end
    nargs = nargs-1;
    varargin(1) = [];
end


if nargs == 0 % there was only a message
    nargs = 1;
    varargin(1) = {0}; % starting number
end

first_msg = false;
ignore = false;
try 
    ignore = pvar.ignore_next;
end

if nargs == 0 || varargin{1} <=0 || isempty(pvar)
    first_msg = true;
    if ~ignore
        pvar.log_level = opt.log_level;
    end
    if pvar.log_level > eidors_msg('log_level');
        return
    end
    if ~ignore
        pvar.final_msg = opt.final_msg;
    end
    if varargin{1} < 0
        pvar.ignore_next = true;
    else
        pvar.ignore_next = false;
    end
    pvar.timerVal = tic;
    pvar.own_time = 0;
    if nargs <= 1
        pvar.nChars = 4;
    elseif varargin{2} == 0
        pvar.nChars = opt.num_digits;
    else
        pvar.numLength = floor(log10(varargin{2})) + 1;
        pvar.nChars = 2 * pvar.numLength + 1;
    end
end

if pvar.log_level > eidors_msg('log_level');
    return
end

if ~isempty(msg) && ~isinf(varargin{1}) && ~(first_msg && ignore);
    fprintf('%s ',msg);
end

if nargs > 0 && isinf(varargin{1})
    if isempty(msg)
        msg = pvar.final_msg;
        if strcmp(msg,'none')
           return
        end
    end
    print_final_msg(msg, pvar.nChars);
    print_time(pvar);
    if eidors_debug('query','progress_msg')
        fprintf('Self-time: %f\n',pvar.own_time);
    end
    return
end
    

if nargs == 1
    percentmsg(varargin{1},first_msg, pvar.nChars);
elseif varargin{2} > 0
    outofmsg(varargin{1},varargin{2},first_msg,...
        pvar.numLength, pvar.nChars);
else
    numbermsg(varargin{1},first_msg, pvar.nChars);
end

pvar.own_time = pvar.own_time + toc(t0);

end


function percentmsg(num, first, N)
    if first && num < 0, return, end
    str = '%3d%%';
    if ~first
        str = [repmat('\b',1,N) str ];
    end
    fprintf(str, round(100*num));
end

function outofmsg(a,b,first,N,T)
    if first && a < 0, return, end
    str = sprintf('%%%dd/%%%dd',N+1,N);
    if ~first
        str = [repmat('\b',1,T+1) str ];
    end
    
    fprintf(str,a,b);
end

function numbermsg(num, first, T)
    if first && num < 0, return, end
    str = sprintf('%%%dd',T);
    if ~first
        str = [repmat('\b',1,T) str];
    end
    fprintf(str, num);
            
end

function print_final_msg(msg, N)
    if isempty(msg) , return, end
    str = [repmat('\b',1,N) '%s'];
    fprintf(str, msg);
end

function print_time(pvar)
    fprintf(' (%.3f s)\n',...
        toc(pvar.timerVal) -  pvar.own_time);
end

function opt = default_opts
    opt.final_msg = 'Done';
    opt.log_level = 2;
    opt.num_digits = 6;
end

function do_unit_test(N)
    eidors_msg('log_level',2);
    eidors_debug on progress_msg
    if nargin == 0
        N = 1:35;
    end
    for n = N
        switch n
            case 1
                progress_msg(0);

            case 2
                progress_msg('Test 2',0);
            case 3
                progress_msg('Test 3');
            case 4
                progress_msg
            case 5
                opt.final_msg = 'YUPPIE!';
                progress_msg(opt);
            case 6
                progress_msg('Test 6', opt);
            case 10
                progress_msg('Test 10',0,10);
            case 11
                progress_msg(0,10);
            case 20
                opt.final_msg = 'Ready !!!';
                progress_msg('Test 20', 0,10,opt);
            case 21
                eidors_msg('log_level',1);
                progress_msg('Test 21', 0,10,opt);
            case 22
                opt.log_level = 1;
                progress_msg('Test 22', 0,10,opt);
            case 23
                opt.log_level = 4;
                eidors_msg('log_level',5);
                progress_msg('Test 23', 0,10,opt);
            case 24
                eidors_msg('log_level',2);
                progress_msg('Test 24', 0,10);
            case 30
                progress_msg('Test 30',0,0)
            otherwise
                continue
        end 
        if n < 10
            for i = 1:10
                pause(.2);
                progress_msg(i/10);
            end
        elseif n >= 30
            for i = 1:10
                pause(.2);
                progress_msg(i^4,0);
            end
        else
            for i = 1:10
                pause(.2);
                progress_msg(i,10);
            end
        end
        if n == 2
           progress_msg('Custom final message', Inf);
        else
           progress_msg(Inf);
        end

    end
    eidors_debug off progress_msg
end
