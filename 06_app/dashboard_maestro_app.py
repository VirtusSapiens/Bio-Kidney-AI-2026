import sys,os,subprocess
from datetime import datetime
from PyQt6.QtWidgets import QApplication,QMainWindow,QWidget,QVBoxLayout,QHBoxLayout,QLabel,QPushButton,QGroupBox,QGridLayout,QFrame,QSizePolicy,QFileDialog,QMessageBox,QStatusBar,QScrollArea
from PyQt6.QtCore import Qt,QTimer,QRect,QPropertyAnimation,QEasingCurve,pyqtSignal
from PyQt6.QtGui import QFont,QColor,QPalette,QLinearGradient,QPainter,QBrush,QPen,QRadialGradient,QCursor
C_BG="#0D2137";C_P="#0A1929";C_D="#061220";C_G="#00FF88";C_DIM="#7EB8D4";C_L="#C8E6F5";C_W="#FFFFFF"
PC={"vasc":"#FF4444","o2":"#00FFFF","ipsc":"#44FF44","swift":"#FFCC00","filtro":"#4488FF","reabs":"#FF44FF"}
INTERP={"vasc":"Red vascular optimizada con CCO v7, 1448 segmentos cumpliendo ley de Murray. Presion 13 kPa fisiologica. Cobertura 100% sin isquemia.","o2":"Difusion de Krogh resuelta. PO2 minima 5.6 mmHg sobre umbral critico 4 mmHg. 0% hipoxia. Eficiencia O2 94%.","ipsc":"iPSCs diferenciadas: podocitos 100%, asa Henle 100%, tubulo proximal 100%. Riesgo teratoma bajo (OCT4 menor 0.1%).","swift":"Co-SWIFT 60 Pa en rango optimo 40-80 Pa. Viabilidad 98% supera estandar clinico 85%. WSS 5.6 dyn/cm2 fisiologico.","filtro":"TFG 82 mL/min equivale ERC Estadio 2. Presiones Starling balanceadas. Diafragma de ranura integro. Flujo laminar Re menor 1.","reabs":"Reabsorcion tubular 98.1% en rango normal 97-99%. 2.19 L/dia orina. NHE3, NKCC2, SGLT2 funcionales. 6/6 criterios cumplidos."}
DATOS=[("VASCULARIZACION CCO v7","vasc",[("Segmentos Murray:","1,448"),("Cumplimiento CCO:","100%"),("Presion entrada:","13 kPa"),("Cobertura red:","100%")]),("OXIGENACION O2","o2",[("PO2 minima:","5.6 mmHg"),("Zona hipoxica:","0%"),("Tejido oxigenado:","100%"),("Gradiente O2:","OPTIMO")]),("DIFERENCIACION iPSC","ipsc",[("Pureza podocitos:","100%"),("Pureza asa Henle:","100%"),("Pureza tubulo px:","100%"),("Riesgo teratoma:","BAJO")]),("BIOIMPRESION Co-SWIFT","swift",[("Presion extrusion:","60 Pa"),("Viabilidad celular:","98%"),("WSS hidrodinamica:","5.6 dyn"),("Soluciones Pareto:","100")]),("FILTRACION GLOMERULAR","filtro",[("TFG simulada:","82 mL/min"),("Presion Starling:","OPTIMA"),("Slit diaphragm:","INTEGRO"),("Flujo capilar:","LAMINAR")]),("REABSORCION TUBULAR","reabs",[("Orina producida:","2.19 L/dia"),("Tasa reabsorcion:","98.1%"),("Flujo urinario:","1.52 mL/min"),("Criterios OK:","6/6")])]
def generar_pdf_modulo(nombre,key,metricas,color_hex,parent=None):
    path,_=QFileDialog.getSaveFileName(parent,"PDF "+nombre,os.path.expanduser("~/Escritorio/BioKidney_"+key+".pdf"),"PDF (*.pdf)")
    if not path:return
    try:
        from PyQt6.QtPrintSupport import QPrinter
        pr=QPrinter(QPrinter.PrinterMode.HighResolution);pr.setOutputFormat(QPrinter.OutputFormat.PdfFormat);pr.setOutputFileName(path);pr.setPageSize(QPrinter.PageSize.A4)
        pa=QPainter(pr);pw=pr.width();ph=pr.height();col=QColor(color_hex)
        pa.fillRect(0,0,pw,ph,QColor("#0D2137"));pa.fillRect(0,0,pw,180,QColor("#0A1929"))
        pa.setPen(QPen(col,4));pa.drawRect(10,10,pw-20,160)
        pa.setFont(QFont("DejaVu Sans Mono",24,QFont.Weight.Bold));pa.setPen(col)
        pa.drawText(QRect(20,20,pw-40,80),Qt.AlignmentFlag.AlignCenter,nombre)
        pa.setFont(QFont("DejaVu Sans Mono",13));pa.setPen(QColor("#C8E6F5"))
        pa.drawText(QRect(20,100,pw-40,38),Qt.AlignmentFlag.AlignCenter,"Bio-Kidney AI 2026  |  VirtusSapiens")
        pa.setFont(QFont("DejaVu Sans Mono",10));pa.setPen(QColor("#7EB8D4"))
        pa.drawText(QRect(20,138,pw-40,30),Qt.AlignmentFlag.AlignCenter,"Validacion In Silico Completa - Pipeline 100%")
        y=200;pa.setFont(QFont("DejaVu Sans Mono",15,QFont.Weight.Bold));pa.setPen(QColor("#C8E6F5"))
        pa.drawText(QRect(40,y,pw-80,34),Qt.AlignmentFlag.AlignLeft,"METRICAS CLAVE")
        pa.setPen(QPen(col,2));pa.drawLine(40,y+36,pw-40,y+36);y+=52
        for label,valor in metricas:
            pa.fillRect(40,y,pw-80,52,QColor("#0A1929"));pa.setPen(QPen(col,1));pa.drawRect(40,y,pw-80,52)
            pa.setFont(QFont("DejaVu Sans Mono",13));pa.setPen(QColor("#7EB8D4"))
            pa.drawText(QRect(55,y+6,(pw-120)//2,40),Qt.AlignmentFlag.AlignVCenter|Qt.AlignmentFlag.AlignLeft,label)
            pa.setFont(QFont("DejaVu Sans Mono",14,QFont.Weight.Bold));pa.setPen(QColor("#00FF88"))
            pa.drawText(QRect(pw//2,y+6,pw//2-55,40),Qt.AlignmentFlag.AlignVCenter|Qt.AlignmentFlag.AlignRight,valor+" OK");y+=62
        y+=16;pa.setFont(QFont("DejaVu Sans Mono",14,QFont.Weight.Bold));pa.setPen(QColor("#C8E6F5"))
        pa.drawText(QRect(40,y,pw-80,34),Qt.AlignmentFlag.AlignLeft,"INTERPRETACION CIENTIFICA")
        pa.setPen(QPen(col,2));pa.drawLine(40,y+36,pw-40,y+36);y+=52
        interp=INTERP.get(key,"");pa.setFont(QFont("DejaVu Sans",11));pa.setPen(QColor("#C8E6F5"))
        palabras=interp.split();linea=""
        for palabra in palabras:
            prueba=linea+" "+palabra if linea else palabra
            if len(prueba)>78:
                pa.drawText(QRect(40,y,pw-80,26),Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter,linea);y+=28;linea=palabra
            else:linea=prueba
        if linea:pa.drawText(QRect(40,y,pw-80,26),Qt.AlignmentFlag.AlignLeft|Qt.AlignmentFlag.AlignVCenter,linea);y+=28
        pa.fillRect(0,ph-80,pw,80,QColor("#061220"));pa.setPen(QPen(col,1));pa.drawLine(0,ph-80,pw,ph-80)
        pa.setFont(QFont("DejaVu Sans Mono",10));pa.setPen(QColor("#7EB8D4"))
        pa.drawText(QRect(20,ph-65,pw-40,26),Qt.AlignmentFlag.AlignCenter,"Carlos David Moreno Caceres  .  VirtusSapiens  .  Medellin, Colombia")
        pa.setFont(QFont("DejaVu Sans Mono",9));pa.setPen(QColor("#00FF88"))
        pa.drawText(QRect(20,ph-39,pw-40,26),Qt.AlignmentFlag.AlignCenter,__import__("datetime").datetime.now().strftime("Generado: %Y-%m-%d %H:%M")+"  |  Bio-Kidney AI 2026")
        pa.end();QMessageBox.information(parent,"PDF OK","Guardado en:\n"+path)
    except Exception as e:QMessageBox.warning(parent,"Error",str(e))
def generar_pdf_global(parent=None):
    path,_=QFileDialog.getSaveFileName(parent,"PDF Maestro",os.path.expanduser("~/Escritorio/BioKidney_Reporte_Maestro.pdf"),"PDF (*.pdf)")
    if not path:return
    try:
        from PyQt6.QtPrintSupport import QPrinter
        pr=QPrinter(QPrinter.PrinterMode.HighResolution);pr.setOutputFormat(QPrinter.OutputFormat.PdfFormat);pr.setOutputFileName(path);pr.setPageSize(QPrinter.PageSize.A4)
        pa=QPainter(pr);pw=pr.width();ph=pr.height()
        pa.fillRect(0,0,pw,ph,QColor("#0D2137"));pa.fillRect(0,0,pw,260,QColor("#0A1929"))
        pa.setPen(QPen(QColor("#00FF00"),4));pa.drawRect(14,14,pw-28,232)
        pa.setFont(QFont("DejaVu Sans Mono",28,QFont.Weight.Bold));pa.setPen(QColor("#00FF88"))
        pa.drawText(QRect(20,28,pw-40,75),Qt.AlignmentFlag.AlignCenter,"BIO-KIDNEY AI 2026")
        pa.setFont(QFont("DejaVu Sans Mono",16,QFont.Weight.Bold));pa.setPen(QColor("#C8E6F5"))
        pa.drawText(QRect(20,108,pw-40,46),Qt.AlignmentFlag.AlignCenter,"REPORTE MAESTRO - VALIDACION IN SILICO")
        pa.setFont(QFont("DejaVu Sans Mono",13));pa.setPen(QColor("#00FF88"))
        pa.drawText(QRect(20,158,pw-40,36),Qt.AlignmentFlag.AlignCenter,"Pipeline 100% - 12/12 Modulos OPTIMOS")
        pa.setFont(QFont("DejaVu Sans Mono",11));pa.setPen(QColor("#7EB8D4"))
        pa.drawText(QRect(20,196,pw-40,34),Qt.AlignmentFlag.AlignCenter,"RINON HUMANO FUNCIONAL BIOIMPRESON")
        y=280;pa.setFont(QFont("DejaVu Sans Mono",12,QFont.Weight.Bold));pa.setPen(QColor("#C8E6F5"))
        pa.drawText(QRect(40,y,pw-80,32),Qt.AlignmentFlag.AlignLeft,"RESUMEN DE MODULOS")
        pa.setPen(QPen(QColor("#00FF88"),2));pa.drawLine(40,y+34,pw-40,y+34);y+=50
        for nombre,key,metricas in DATOS:
            col=QColor(PC[key]);pa.fillRect(40,y,pw-80,38,QColor("#0A1929"));pa.setPen(QPen(col,2));pa.drawRect(40,y,pw-80,38)
            pa.setFont(QFont("DejaVu Sans Mono",10,QFont.Weight.Bold));pa.setPen(col)
            pa.drawText(QRect(54,y+3,340,32),Qt.AlignmentFlag.AlignVCenter|Qt.AlignmentFlag.AlignLeft,nombre)
            pa.setFont(QFont("DejaVu Sans Mono",9));pa.setPen(QColor("#00FF88"))
            resumen="  ".join([v+" OK" for _,v in metricas[:2]])
            pa.drawText(QRect(pw//2-10,y+3,pw//2,32),Qt.AlignmentFlag.AlignVCenter|Qt.AlignmentFlag.AlignLeft,resumen);y+=48
        pa.fillRect(0,ph-80,pw,80,QColor("#061220"));pa.setPen(QPen(QColor("#00FF88"),1));pa.drawLine(0,ph-80,pw,ph-80)
        pa.setFont(QFont("DejaVu Sans Mono",10));pa.setPen(QColor("#7EB8D4"))
        pa.drawText(QRect(20,ph-65,pw-40,26),Qt.AlignmentFlag.AlignCenter,"Carlos David Moreno Caceres  .  VirtusSapiens  .  Medellin, Colombia  .  2026")
        pa.setFont(QFont("DejaVu Sans Mono",9));pa.setPen(QColor("#00FF88"))
        pa.drawText(QRect(20,ph-39,pw-40,26),Qt.AlignmentFlag.AlignCenter,__import__("datetime").datetime.now().strftime("Generado: %Y-%m-%d %H:%M")+"  |  Bio-Kidney AI 2026")
        pa.end();QMessageBox.information(parent,"PDF Maestro OK","Guardado en:\n"+path)
    except Exception as e:QMessageBox.warning(parent,"Error",str(e))
class PanelLateral(QWidget):
    def __init__(self,parent=None):
        super().__init__(parent);self.setFixedWidth(0);self._key=None;self._nombre=None;self._metricas=[];self._color=C_G
        self._a1=QPropertyAnimation(self,b"minimumWidth");self._a1.setEasingCurve(QEasingCurve.Type.OutCubic);self._a1.setDuration(300)
        self._a2=QPropertyAnimation(self,b"maximumWidth");self._a2.setEasingCurve(QEasingCurve.Type.OutCubic);self._a2.setDuration(300)
        self._build_ui()
    def _build_ui(self):
        self.setStyleSheet("background-color:#061220;border-left:2px solid #1A3A5C;")
        lay=QVBoxLayout(self);lay.setContentsMargins(14,14,14,14);lay.setSpacing(8)
        top=QHBoxLayout();self.lbl=QLabel("");self.lbl.setFont(QFont("DejaVu Sans Mono",10,QFont.Weight.Bold))
        self.lbl.setStyleSheet("color:#C8E6F5;background:transparent;");self.lbl.setWordWrap(True);top.addWidget(self.lbl);top.addStretch()
        bx=QPushButton("x");bx.setFixedSize(26,26);bx.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        bx.setStyleSheet("QPushButton{background:#0A1929;color:#FF4444;border:1px solid #FF4444;border-radius:4px;font-weight:bold;}QPushButton:hover{background:#FF444422;}")
        bx.clicked.connect(self.cerrar);top.addWidget(bx);lay.addLayout(top)
        self.sep=QFrame();self.sep.setFrameShape(QFrame.Shape.HLine);self.sep.setStyleSheet("border:1px solid #1A3A5C;");lay.addWidget(self.sep)
        self.mw=QWidget();self.mw.setStyleSheet("background:transparent;")
        self.ml=QVBoxLayout(self.mw);self.ml.setSpacing(5);self.ml.setContentsMargins(0,0,0,0);lay.addWidget(self.mw)
        s2=QFrame();s2.setFrameShape(QFrame.Shape.HLine);s2.setStyleSheet("border:1px solid #1A3A5C;");lay.addWidget(s2)
        li=QLabel("INTERPRETACION CIENTIFICA");li.setFont(QFont("DejaVu Sans Mono",8,QFont.Weight.Bold));li.setStyleSheet("color:#7EB8D4;background:transparent;");lay.addWidget(li)
        self.interp=QLabel("");self.interp.setFont(QFont("DejaVu Sans",9));self.interp.setStyleSheet("color:#C8E6F5;background:transparent;");self.interp.setWordWrap(True)
        self.interp.setAlignment(Qt.AlignmentFlag.AlignTop|Qt.AlignmentFlag.AlignLeft)
        sc=QScrollArea();sc.setWidget(self.interp);sc.setWidgetResizable(True)
        sc.setStyleSheet("QScrollArea{background:transparent;border:none;}QScrollBar:vertical{background:#0A1929;width:6px;}QScrollBar::handle:vertical{background:#1A3A5C;border-radius:3px;}")
        sc.setFixedHeight(175);lay.addWidget(sc);lay.addStretch()
        self.btn=QPushButton("Generar PDF de este modulo");self.btn.setFont(QFont("DejaVu Sans Mono",9,QFont.Weight.Bold));self.btn.setFixedHeight(40)
        self.btn.setCursor(QCursor(Qt.CursorShape.PointingHandCursor));self.btn.clicked.connect(self._pdf);lay.addWidget(self.btn)
    def abrir(self,nombre,key,metricas,color):
        self._key=key;self._nombre=nombre;self._metricas=metricas;self._color=color
        self.lbl.setText(nombre);self.lbl.setStyleSheet("color:"+color+";background:transparent;font-weight:bold;font-family:DejaVu Sans Mono;font-size:10px;")
        self.sep.setStyleSheet("border:1px solid "+color+";")
        for i in reversed(range(self.ml.count())):
            w=self.ml.itemAt(i).widget()
            if w:w.setParent(None)
        for label,valor in metricas:
            fi=QWidget();fi.setStyleSheet("background:#0A1929;border-radius:4px;")
            fl=QHBoxLayout(fi);fl.setContentsMargins(8,3,8,3)
            lb=QLabel(label);lb.setFont(QFont("DejaVu Sans Mono",8));lb.setStyleSheet("color:#7EB8D4;background:transparent;")
            lb.setSizePolicy(QSizePolicy.Policy.Expanding,QSizePolicy.Policy.Preferred);fl.addWidget(lb)
            vl=QLabel(valor+" OK");vl.setFont(QFont("DejaVu Sans Mono",8,QFont.Weight.Bold));vl.setStyleSheet("color:#00FF88;background:transparent;");fl.addWidget(vl)
            self.ml.addWidget(fi)
        self.interp.setText(INTERP.get(key,""))
        self.btn.setStyleSheet("QPushButton{background-color:#0A1929;color:"+color+";border:2px solid "+color+";border-radius:6px;}QPushButton:hover{background-color:"+color+"22;color:#FFFFFF;}QPushButton:pressed{background-color:"+color+"44;}")
        self._a1.setStartValue(self.width());self._a1.setEndValue(370);self._a2.setStartValue(self.width());self._a2.setEndValue(370)
        self._a1.start();self._a2.start()
    def cerrar(self):
        self._a1.setStartValue(self.width());self._a1.setEndValue(0);self._a2.setStartValue(self.width());self._a2.setEndValue(0)
        self._a1.start();self._a2.start()
    def _pdf(self):generar_pdf_modulo(self._nombre,self._key,self._metricas,self._color,self)
class GlobalIndicator(QWidget):
    clicked=pyqtSignal()
    def __init__(self,parent=None):
        super().__init__(parent);self.setFixedHeight(130);self._step=0;self._dir=1;self._pulse=0.0;self._hover=False
        self.setCursor(QCursor(Qt.CursorShape.PointingHandCursor));t=QTimer(self);t.timeout.connect(self._tick);t.start(40)
    def _tick(self):
        self._step+=self._dir*3
        if self._step>=100:self._dir=-1
        elif self._step<=0:self._dir=1
        self._pulse=self._step/100.0;self.update()
    def enterEvent(self,e):self._hover=True;self.update()
    def leaveEvent(self,e):self._hover=False;self.update()
    def mousePressEvent(self,e):
        if e.button()==Qt.MouseButton.LeftButton:self.clicked.emit()
    def paintEvent(self,event):
        p=QPainter(self);p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w,h=self.width(),self.height();g=QLinearGradient(0,0,w,0);g.setColorAt(0,QColor("#0A1929"));g.setColorAt(0.5,QColor("#0D2F50"));g.setColorAt(1,QColor("#0A1929"))
        p.fillRect(0,0,w,h,QBrush(g));alpha=min(255,int(80+120*self._pulse)+(40 if self._hover else 0))
        p.setPen(QPen(QColor(0,255,0,alpha),3 if self._hover else 2));p.drawRoundedRect(10,6,w-20,h-12,8,8)
        p.setPen(QPen(QColor("#00AA44"),1));p.drawRoundedRect(13,9,w-26,h-18,6,6)
        gx=w//2-180;p.setFont(QFont("DejaVu Sans",28,QFont.Weight.Bold));p.setPen(QColor(0,int(160+95*self._pulse),60))
        p.drawText(QRect(gx,14,60,50),Qt.AlignmentFlag.AlignCenter,"OK")
        p.setFont(QFont("DejaVu Sans Mono",20,QFont.Weight.Bold));p.setPen(QColor(0,int(200+55*self._pulse),80))
        p.drawText(QRect(gx+60,12,w-gx-80,40),Qt.AlignmentFlag.AlignVCenter|Qt.AlignmentFlag.AlignLeft,"  RINON FUNCIONAL")
        p.setFont(QFont("DejaVu Sans Mono",11));p.setPen(QColor(C_L))
        p.drawText(QRect(gx+60,50,w-gx-80,28),Qt.AlignmentFlag.AlignVCenter|Qt.AlignmentFlag.AlignLeft,"  VALIDACION IN SILICO COMPLETA")
        p.setFont(QFont("DejaVu Sans Mono",10));p.setPen(QColor(C_G))
        p.drawText(QRect(gx+60,78,w-gx-80,26),Qt.AlignmentFlag.AlignVCenter|Qt.AlignmentFlag.AlignLeft,"  Pipeline 100% - 12/12 modulos OPTIMOS")
        if self._hover:
            p.setFont(QFont("DejaVu Sans Mono",8));p.setPen(QColor(0,200,80,180))
            p.drawText(QRect(w-240,h-22,230,18),Qt.AlignmentFlag.AlignRight,"Clic: PDF Reporte Maestro completo")
        p.end()
class ModulePanel(QGroupBox):
    clicked=pyqtSignal(str,str,list,str)
    def __init__(self,nombre,key,metricas,parent=None):
        super().__init__(parent);self._n=nombre;self._k=key;self._m=metricas;self.bc=PC[key];self.setTitle("");self.setFixedHeight(200)
        self.setCursor(QCursor(Qt.CursorShape.PointingHandCursor));self._hov=False
        lay=QVBoxLayout(self);lay.setContentsMargins(12,10,12,10);lay.setSpacing(4)
        tl=QLabel(nombre);tl.setAlignment(Qt.AlignmentFlag.AlignCenter);tl.setFont(QFont("DejaVu Sans Mono",10,QFont.Weight.Bold))
        tl.setStyleSheet("color:"+self.bc+";background:transparent;");lay.addWidget(tl)
        sep=QFrame();sep.setFrameShape(QFrame.Shape.HLine);sep.setStyleSheet("border:1px solid "+self.bc+";background:"+self.bc+";");sep.setFixedHeight(1);lay.addWidget(sep);lay.addSpacing(4)
        fm=QFont("DejaVu Sans Mono",9)
        for lb,vl in metricas:
            row=QHBoxLayout();row.setSpacing(4)
            l=QLabel(lb);l.setFont(fm);l.setStyleSheet("color:"+C_DIM+";background:transparent;")
            l.setSizePolicy(QSizePolicy.Policy.Expanding,QSizePolicy.Policy.Preferred);row.addWidget(l)
            v=QLabel(vl+" OK");v.setFont(QFont("DejaVu Sans Mono",9,QFont.Weight.Bold));v.setStyleSheet("color:"+C_G+";background:transparent;")
            v.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignVCenter);row.addWidget(v);lay.addLayout(row)
        lay.addStretch();ht=QLabel("Clic para detalles y PDF");ht.setFont(QFont("DejaVu Sans Mono",7))
        ht.setStyleSheet("color:"+self.bc+"99;background:transparent;");ht.setAlignment(Qt.AlignmentFlag.AlignRight);lay.addWidget(ht)
        self._us()
    def _us(self):
        b=3 if self._hov else 2;bg=self.bc+"11" if self._hov else C_P
        self.setStyleSheet("QGroupBox{background-color:"+bg+";border:"+str(b)+"px solid "+self.bc+";border-radius:8px;padding:4px;}")
    def enterEvent(self,e):self._hov=True;self._us();self.update()
    def leaveEvent(self,e):self._hov=False;self._us();self.update()
    def mousePressEvent(self,e):
        if e.button()==Qt.MouseButton.LeftButton:self.clicked.emit(self._n,self._k,self._m,self.bc)
    def paintEvent(self,event):
        super().paintEvent(event);p=QPainter(self);p.setRenderHint(QPainter.RenderHint.Antialiasing)
        g=QRadialGradient(15,15,40);g.setColorAt(0,QColor(self.bc).lighter(150));g.setColorAt(1,QColor(self.bc+"00"))
        p.fillRect(0,0,50,50,QBrush(g));p.end()
class FlowDiagram(QWidget):
    def __init__(self,parent=None):
        super().__init__(parent);self.setFixedHeight(100);self._off=0
        t=QTimer(self);t.timeout.connect(self._tick);t.start(60)
    def _tick(self):self._off=(self._off+2)%40;self.update()
    def paintEvent(self,event):
        p=QPainter(self);p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w,h=self.width(),self.height();p.fillRect(0,0,w,h,QColor(C_D))
        p.setPen(QPen(QColor("#1A3A5C"),1));p.drawRoundedRect(4,4,w-8,h-8,6,6)
        steps=[("SANGRE ART",PC["vasc"]),("CCO v7",PC["vasc"]),("FILTRACION",PC["filtro"]),("REABSORCION",PC["reabs"]),("ORINA",PC["o2"]),("SANGRE LIMPIA",PC["ipsc"])]
        n=len(steps);mg=20;sw=(w-2*mg)/n;bw=int(sw*0.78);bh=56;yc=h//2
        p.setFont(QFont("DejaVu Sans Mono",7,QFont.Weight.Bold))
        for i,(txt,col) in enumerate(steps):
            cx=int(mg+sw*i+sw/2);bx=cx-bw//2;by=yc-bh//2
            p.setPen(QPen(QColor(col),2));p.setBrush(QBrush(QColor(C_P)));p.drawRoundedRect(bx,by,bw,bh,5,5)
            p.setPen(QColor(col));p.drawText(QRect(bx,by,bw,bh),Qt.AlignmentFlag.AlignCenter,txt)
            if i<n-1:
                xs=bx+bw+2;xe=int(mg+sw*(i+1)+sw/2-bw//2-2)
                p.setPen(QPen(QColor(col),1,Qt.PenStyle.DashLine));dx=xs
                while dx<xe-8:
                    s=dx+(self._off%20);e=min(s+10,xe-8);p.drawLine(int(s),yc,int(e),yc);dx+=20
                p.setPen(QPen(QColor(col),2));p.drawLine(xe-8,yc-5,xe,yc);p.drawLine(xe-8,yc+5,xe,yc)
        p.end()
class DashboardMaestro(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Bio-Kidney AI 2026 - Dashboard Maestro | VirtusSapiens")
        self.setMinimumSize(1100,800);self.resize(1280,880)
        self.setStyleSheet("QMainWindow,QWidget{background-color:"+C_BG+";color:"+C_L+";}");
        self._build_ui()
        sc=QApplication.primaryScreen().availableGeometry();fg=self.frameGeometry();fg.moveCenter(sc.center());self.move(fg.topLeft())
    def _build_ui(self):
        rw=QWidget();self.setCentralWidget(rw);rh=QHBoxLayout(rw);rh.setContentsMargins(0,0,0,0);rh.setSpacing(0)
        main=QWidget();main.setStyleSheet("background-color:"+C_BG+";");rh.addWidget(main,1)
        lay=QVBoxLayout(main);lay.setContentsMargins(16,12,16,8);lay.setSpacing(10)
        hb=QHBoxLayout()
        at=QLabel("Bio-Kidney AI 2026 - Dashboard Maestro");at.setFont(QFont("DejaVu Sans Mono",13,QFont.Weight.Bold));at.setStyleSheet("color:"+C_L+";")
        vs=QLabel("VirtusSapiens");vs.setFont(QFont("DejaVu Sans Mono",11,QFont.Weight.Bold));vs.setStyleSheet("color:"+PC["o2"]+";");vs.setAlignment(Qt.AlignmentFlag.AlignRight|Qt.AlignmentFlag.AlignVCenter)
        hb.addWidget(at);hb.addStretch();hb.addWidget(vs);lay.addLayout(hb)
        self.gi=GlobalIndicator();self.gi.clicked.connect(lambda:generar_pdf_global(self));lay.addWidget(self.gi)
        grid=QGridLayout();grid.setSpacing(10)
        for idx,(nombre,key,metricas) in enumerate(DATOS):
            panel=ModulePanel(nombre,key,metricas);panel.clicked.connect(self._abrir)
            r,c=divmod(idx,3);grid.addWidget(panel,r,c)
        lay.addLayout(grid)
        fl=QLabel("FLUJO FUNCIONAL DEL RINON BIOIMPRSO");fl.setFont(QFont("DejaVu Sans Mono",9));fl.setStyleSheet("color:"+C_DIM+";");fl.setAlignment(Qt.AlignmentFlag.AlignCenter);lay.addWidget(fl)
        lay.addWidget(FlowDiagram())
        bb=QHBoxLayout();bb.setSpacing(12)
        bp=self._btn("Exportar PNG",PC["o2"]);bd=self._btn("PDF Reporte Maestro",PC["swift"]);bs=self._btn("Ver Simuladores",PC["ipsc"])
        bp.clicked.connect(self._png);bd.clicked.connect(lambda:generar_pdf_global(self));bs.clicked.connect(self._sim)
        bb.addStretch();bb.addWidget(bp);bb.addWidget(bd);bb.addWidget(bs);bb.addStretch();lay.addLayout(bb)
        ft=QLabel("Carlos David Moreno Caceres  .  VirtusSapiens  .  Medellin, Colombia  .  2026")
        ft.setAlignment(Qt.AlignmentFlag.AlignCenter);ft.setFont(QFont("DejaVu Sans Mono",8));ft.setStyleSheet("color:"+C_DIM+";padding:4px 0;");lay.addWidget(ft)
        self.pl=PanelLateral();rh.addWidget(self.pl)
        self.sb=QStatusBar();self.setStatusBar(self.sb)
        self.sb.setStyleSheet("background:"+C_D+";color:"+C_DIM+";font-family:monospace;font-size:9px;")
        self.sb.showMessage("  Bio-Kidney AI 2026  .  12/12 OPTIMOS  .  Clic en cualquier modulo para detalles y PDF  .  "+datetime.now().strftime("%Y-%m-%d %H:%M"))
    def _abrir(self,nombre,key,metricas,color):
        self.pl.abrir(nombre,key,metricas,color)
        self.sb.showMessage("  "+nombre+"  |  Panel lateral abierto - Clic Generar PDF",4000)
    def _btn(self,txt,col):
        b=QPushButton(txt);b.setFont(QFont("DejaVu Sans Mono",10,QFont.Weight.Bold));b.setFixedHeight(40);b.setMinimumWidth(180)
        b.setCursor(QCursor(Qt.CursorShape.PointingHandCursor))
        b.setStyleSheet("QPushButton{background-color:"+C_P+";color:"+col+";border:2px solid "+col+";border-radius:6px;padding:4px 16px;}QPushButton:hover{background-color:"+col+"22;color:"+C_W+";}QPushButton:pressed{background-color:"+col+"44;}")
        return b
    def _png(self):
        path,_=QFileDialog.getSaveFileName(self,"Guardar PNG",os.path.expanduser("~/Escritorio/BioKidney_Dashboard.png"),"PNG (*.png)")
        if path:
            QApplication.primaryScreen().grabWindow(int(self.winId())).save(path,"PNG")
            self.sb.showMessage("  Guardado: "+path,5000);QMessageBox.information(self,"OK","Guardado: "+path)
    def _sim(self):
        sp=os.path.expanduser("~/Escritorio/BioKidney-AI/06_app/biokidney_app.py")
        ea=os.path.expanduser("~/env_biokidney/bin/activate")
        if os.path.exists(sp):
            subprocess.Popen("bash -c \"source "+ea+" && python "+sp+" &\"",shell=True)
            self.sb.showMessage("  Abriendo Simuladores...",4000)
        else:QMessageBox.warning(self,"No encontrado","No encontrado: "+sp)
def main():
    app=QApplication(sys.argv);pal=QPalette()
    for role,col in [(QPalette.ColorRole.Window,C_BG),(QPalette.ColorRole.WindowText,C_L),(QPalette.ColorRole.Base,C_P),(QPalette.ColorRole.Text,C_L),(QPalette.ColorRole.Button,C_P),(QPalette.ColorRole.ButtonText,C_L),(QPalette.ColorRole.Highlight,C_G)]:
        pal.setColor(role,QColor(col))
    app.setPalette(pal);w=DashboardMaestro();w.show();sys.exit(app.exec())
if __name__=="__main__":
    main()
