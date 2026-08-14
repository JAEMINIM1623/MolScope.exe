@echo off
chcp 65001 >nul
REM ============================================================
REM  MolScope Windows 실행파일 빌드 스크립트
REM  요구사항: Python 3.10+ 이 설치된 Windows PC
REM  사용법: 이 파일을 더블클릭 (또는 cmd 에서 실행)
REM  결과물: dist\MolScope.exe
REM ============================================================
cd /d "%~dp0"

echo [1/3] 의존성 설치...
python -m pip install --upgrade pip
python -m pip install -r requirements.txt pyinstaller
if errorlevel 1 goto :fail

echo [2/3] 코어 동작 확인...
python -c "from molscope import core; ok, ver = core.rdkit_available(); assert ok, ver; print('RDKit', ver, '| ChemDraw .cdx:', core.chemdraw_cdx_supported())"
if errorlevel 1 goto :fail

echo [3/3] EXE 빌드 (수 분 소요)...
pyinstaller --noconfirm --clean --onefile --windowed ^
  --name MolScope ^
  --collect-all rdkit ^
  --collect-submodules molscope ^
  --hidden-import PIL._tkinter_finder ^
  MolScope.py
if errorlevel 1 goto :fail

echo.
echo  빌드 완료:  dist\MolScope.exe
echo  이 파일 하나만 복사하면 다른 PC 에서도 실행됩니다 (Python 불필요).
pause
exit /b 0

:fail
echo.
echo  빌드 실패. 위 오류 메시지를 확인하세요.
pause
exit /b 1
