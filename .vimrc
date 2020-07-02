" Enable syntax
:syntax on

" set clipboard=unnamed,unnamedplus
set hlsearch
" map :noh, which clear search hilight
nnoremap <silent> ,<space> :nohlsearch<CR>

" shortcut for leaving insert mode
inoremap jj <ESC>

" show line numbers and lenght
set number
" set tw=79
set colorcolumn=100

" status line
" display file name
set laststatus=2
" set statusline=%F "tail of the filename
set statusline=
set statusline+=%#DiffAdd#%{(mode()=='n')?'\ \ NORMAL\ ':''}
set statusline+=%#DiffChange#%{(mode()=='i')?'\ \ INSERT\ ':''}
set statusline+=%#DiffDelete#%{(mode()=='r')?'\ \ RPLACE\ ':''}
set statusline+=%#Cursor#%{(mode()=='v')?'\ \ VISUAL\ ':''}
set statusline+=\ %n\           " buffer number
set statusline+=%#Visual#       " colour
set statusline+=%{&paste?'\ PASTE\ ':''}
set statusline+=%{&spell?'\ SPELL\ ':''}
set statusline+=%#CursorIM#     " colour
set statusline+=%R                        " readonly flag
set statusline+=%M                        " modified [+] flag
set statusline+=%#Cursor#               " colour
set statusline+=%#CursorLine#     " colour
set statusline+=\ %F\                   " long file name isntead of %t for short
set statusline+=%=                          " right align
set statusline+=%#CursorLine#   " colour
set statusline+=\ %Y\                   " file type
set statusline+=%#CursorIM#     " colour
set statusline+=\ %3l:%-2c\         " line + column
set statusline+=%#Cursor#       " colour
set statusline+=\ %3p%%\                " percentage

" following PEP8 indentation 4 python
autocmd BufNewFile,BufRead *.py
\ set tabstop=4
\ | set softtabstop=4
\ | set shiftwidth=4
\ | set shiftround
\ | set expandtab
\ | set fileformat=unix

" folding codes
set foldmethod=indent
set foldlevel=99


" Setup Pathogen to manage your plugins
" mkdir -p ~/.vim/autoload ~/.vim/bundle
" curl -so ~/.vim/autoload/pathogen.vim
" https://raw.githubusercontent.com/tpope/vim-pathogen/master/autoload/pathogen.vim
" Now you can install any plugin into a .vim/bundle/plugin-name/ folder
execute pathogen#infect()

" installing and set up syntastic following: https://github.com/vim-syntastic/syntastic
set statusline+=%#warningmsg#
set statusline+=%{SyntasticStatuslineFlag()}
set statusline+=%*

let g:syntastic_always_populate_loc_list = 1
let g:syntastic_auto_loc_list = 1
let g:syntastic_check_on_open = 1
let g:syntastic_check_on_wq = 0
" let g:syntastic_python_python_exec = 'python3'
" let g:syntastic_python_checkers = ['python']

" Show whitespace
" " MUST be inserted BEFORE the colorscheme command
autocmd ColorScheme * highlight ExtraWhitespace ctermbg=red guibg=red
au InsertLeave * match ExtraWhitespace /\s\+$/

" Color scheme
" wombat:
" mkdir -p ~/.vim/colors && cd ~/.vim/colors
" wget -O wombat256mod.vim http://www.vim.org/scripts/download_script.php?src_id=13400
" set t_Co=256

" gruvbox:
" installing gruvbox colorscheme: following:https://github.com/morhetz/gruvbox/wiki/Installation
color gruvbox
:set bg=dark


"" fixing backspace in new install of vim
"set backspace=indent,eol,start


" " PLUMED STUFFS
" " This allows including the proper PLUMED syntax file:
" ":let &runtimepath.=','.$PLUMED_VIMPATH
" " The former command requires PLUMED_VIMPATH to be set. Alternatively, use
" this:
:let &runtimepath.=',/usr/users/cmb/shared/opt/plumed/v2.5.1/lib/plumed/vim'
" " properly adjusted to the path where PLUMED is installed.
" " This makes autocompletion work in the expected way:
:set completeopt=longest,menuone
" " This enables bindings of F2/F3/F4 to plumed specific commands:
:let plumed_shortcuts=1

" TAB to change the hilighted column
: map <F3> :PMinus<CR>
: map <F4> :PPlus<CR>

" 4smaug
let hostname=system("hostname -s | tr -d '\n'")
if hostname == 'smaug'
    set directory=/home/users/jeremy/.vim/swapfiles
    set backupdir=/home/users/jeremy/.vim/tmp
endif

" system clipboard for wayland display server, wl-clipboard must be installed first: sudo apt-get update, sudo apt-get install wl-clipboard 
let serverdisp = $XDG_SESSION_TYPE
if serverdisp == 'wayland'
    xnoremap "+y y:call system("wl-copy", @")<cr>
    nnoremap "+p :let @"=substitute(system("wl-paste --no-newline"), '<C-v><C-m>', '', 'g')<cr>p
    nnoremap "*p :let @"=substitute(system("wl-paste --no-newline --primary"), '<C-v><C-m>', '', 'g')<cr>p
endif
