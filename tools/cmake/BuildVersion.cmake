#building version number in variable _versionString

MACRO (GETDATE TODAY)
    IF (WIN32)
        EXECUTE_PROCESS(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE ${TODAY})
        string(REGEX REPLACE "(..)/(..)/(....).*" "\\3\\2\\1" ${TODAY} ${${TODAY}})
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%Y-%m-%d" OUTPUT_VARIABLE ${TODAY})
        string(REGEX REPLACE "(....)-(..)-(..).*" "\\1\\2\\3" ${TODAY} ${${TODAY}})
    ELSE (WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${TODAY} 000000)
    ENDIF (WIN32)
ENDMACRO (GETDATE)

MACRO(BuildVersion)
SET(_versionString "${VERSION_MAJOR}.${VERSION_MINOR}${VERSION_PATCH}")
ENDMACRO(BuildVersion)
