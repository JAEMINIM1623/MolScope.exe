# -*- coding: utf-8 -*-
"""MolScope 실행 진입점. PyInstaller 는 이 파일을 빌드 대상으로 삼는다."""
import multiprocessing

if __name__ == "__main__":
    multiprocessing.freeze_support()   # PyInstaller onefile 필수
    from molscope.gui import main
    main()
