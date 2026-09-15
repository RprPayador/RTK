#include "RTK_Structs.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// 读取配置文件函数 (简明 C 风格，适合课程项目)
bool ReadRTKConfigInfo(const char FName[], ROVERCFGINFO& cfg) {
    FILE* fp = fopen(FName, "r");
    if (!fp) {
        printf("Open config file '%s' failed.\n", FName);
        return false;
    }

    char line[256];
    char key[64], val[128];

    while (fgets(line, sizeof(line), fp)) {
        // 去除行末注释和换行符
        char* pComment = strpbrk(line, "#;\r\n");
        if (pComment) *pComment = '\0';

        // 查找等号，没有等号则跳过（如空行、纯注释行或 [Section] 分节行）
        char* pEq = strchr(line, '=');
        if (!pEq) continue;

        // 将等号替换为空格，利用 sscanf 快速提取 key 和 val
        *pEq = ' ';
        if (sscanf(line, "%s %s", key, val) < 2) continue;

        // 根据 key 匹配对应结构体成员
        if (strcmp(key, "IsRTKMode") == 0)          cfg.IsRTKMode = (short)atoi(val);
        else if (strcmp(key, "IsFileData") == 0)         cfg.IsFileData = (short)atoi(val);
        else if (strcmp(key, "RTKProcMode") == 0)   cfg.RTKProcMode = (short)atoi(val);
        else if (strcmp(key, "RovPort") == 0)       cfg.RovPort = atoi(val);
        else if (strcmp(key, "RovBaud") == 0)       cfg.RovBaud = atoi(val);
        else if (strcmp(key, "BasNetIP") == 0)      strcpy(cfg.BasNetIP, val);
        else if (strcmp(key, "RovNetIP") == 0)      strcpy(cfg.RovNetIP, val);
        else if (strcmp(key, "BasNetPort") == 0)    cfg.BasNetPort = (short)atoi(val);
        else if (strcmp(key, "RovNetPort") == 0)    cfg.RovNetPort = (short)atoi(val);
        else if (strcmp(key, "CodeNoise") == 0)     cfg.CodeNoise = atof(val);
        else if (strcmp(key, "CPNoise") == 0)       cfg.CPNoise = atof(val);
        else if (strcmp(key, "ElevThreshold") == 0) cfg.ElevThreshold = atof(val);
        else if (strcmp(key, "RatioThres") == 0)    cfg.RatioThres = atof(val);
        else if (strcmp(key, "BasObsDatFile") == 0) strcpy(cfg.BasObsDatFile, val);
        else if (strcmp(key, "RovObsDatFile") == 0) strcpy(cfg.RovObsDatFile, val);
        else if (strcmp(key, "ResFile") == 0)       strcpy(cfg.ResFile, val);
    }

    fclose(fp);
    printf("Config file '%s' loaded successfully.\n", FName);
    return true;
}
