@echo off
setlocal enabledelayedexpansion

:: Name der Ausgabedatei
set OUTPUT=SHGO_current_state.txt

:: Überschrift
echo === SHGO.jl - Aktueller Code-Stand === > "%OUTPUT%"
echo Generiert am %date% um %time% >> "%OUTPUT%"
echo. >> "%OUTPUT%"
echo ------------------------Dateitrennzeichen--------------------------------------------- >> "%OUTPUT%"

:: Alle relevanten Dateien rekursiv verarbeiten
for /r %%F in (*.jl, Project.toml, README.md, *.yml) do (
    if exist "%%F" (
        echo. >> "%OUTPUT%"
        echo "Datei: %%F" >> "%OUTPUT%"
        echo. >> "%OUTPUT%"
        type "%%F" >> "%OUTPUT%"
        echo. >> "%OUTPUT%"
        echo ------------------------Dateitrennzeichen--------------------------------------------- >> "%OUTPUT%"
    )
)

echo Fertig! Alle Dateien wurden in "%OUTPUT%" zusammengefasst.
echo Du kannst diese Datei jetzt kopieren und in ein neues Chat-Fenster einfügen.