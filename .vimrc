" Enable syntax
:syntax on
" set clipboard^=unnamed,unnamedplus

" disable last hilighting
:noh

" shortcut for leaving insert mode
inoremap jj <ESC>

" show line numbers and lenght
set number
" set tw=79
set colorcolumn=80

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
let g:syntastic_python_python_exec = 'python3'
let g:syntastic_python_checkers = ['python']

" Show whitespace
" " MUST be inserted BEFORE the colorscheme command
autocmd ColorScheme * highlight ExtraWhitespace ctermbg=red guibg=red
au InsertLeave * match ExtraWhitespace /\s\+$/

" Color scheme
" mkdir -p ~/.vim/colors && cd ~/.vim/colors
" wget -O wombat256mod.vim http://www.vim.org/scripts/download_script.php?src_id=13400
set t_Co=256


" installing gruvbox colorscheme: following:https://github.com/morhetz/gruvbox/wiki/Installation
color gruvbox
:set bg=dark

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
