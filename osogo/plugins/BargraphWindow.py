#!/usr/bin/env python2

import string

from OsogoPluginWindow import *
from ecell.ecssupport import *

import Numeric
import GTK
import operator

class BargraphWindow( OsogoPluginWindow ):
    
	def __init__( self, dirname, data, pluginmanager, root=None ):

		OsogoPluginWindow.__init__( self, dirname, data, pluginmanager, root )
		self.theSession = pluginmanager.theSession
		aFullPNString = createFullPNString( self.theFullPN() )
		aValue = self.theSession.theSimulator.getProperty( aFullPNString )

		#if operator.isNumberType( aValue[0] ):
		if operator.isNumberType( aValue ):
			self.thePluginManager.appendInstance( self )   
			#self.initialize()
			# -------------------------------------------------
			self['toolbar5'].set_style( GTK.TOOLBAR_ICONS )
			self['toolbar6'].set_style( GTK.TOOLBAR_ICONS )
			self['toolbar5'].set_button_relief( GTK.RELIEF_HALF )
			self['toolbar6'].set_button_relief( GTK.RELIEF_HALF )        
        
			self.pull = 0
			self.thePositiveFlag = 1
			self.theAutoChangeFlag = 1
			self.theActualValue = 0
			self.theBarLength = 0
			self.theMultiplier = 0
        
			self.addHandlers( { \
		                   'on_add_button_clicked'      : self.updateByAddbutton,
		                   'on_subtract_button_clicked' : self.updateBySubtractbutton,
		                   'multiplier_entry_activate'  : self.updateByTextentry,
		                   'auto_button_toggled'        : self.updateByAutoButton ,
		                   'window_exit'                : self.exit })
        
			self.theIDEntry = self.getWidget( "property_id_label" )
			self.theMultiplier1Entry = self.getWidget("multiplier1_label")
			self.update()

			# -------------------------------------------------

		else:
			self.theSession.printMessage( "%s: not numerical data\n" % aFullPNString )

		if len( self.theFullPNList() ) > 1:
			self.addPopupMenu(1,1,1)
		else:
			self.addPopupMenu(0,1,1)


	#def initialize( self ):
    
	#	self['toolbar5'].set_style( GTK.TOOLBAR_ICONS )
	#	self['toolbar6'].set_style( GTK.TOOLBAR_ICONS )
	#	self['toolbar5'].set_button_relief( GTK.RELIEF_HALF )
	#	self['toolbar6'].set_button_relief( GTK.RELIEF_HALF )        
        
	#	self.pull = 0
	#	self.thePositiveFlag = 1
	#	self.theAutoChangeFlag = 1
	#	self.theActualValue = 0
	#	self.theBarLength = 0
	#	self.theMultiplier = 0
        
	#	self.addHandlers( { \
	#	                   'on_add_button_clicked'      : self.updateByAddbutton,
	#	                   'on_subtract_button_clicked' : self.updateBySubtractbutton,
	#	                   'multiplier_entry_activate'  : self.updateByTextentry,
	#	                   'auto_button_toggled'        : self.updateByAutoButton ,
	#	                   'window_exit'                : self.exit })
        
	#	self.theIDEntry = self.getWidget( "property_id_label" )
	#	self.theMultiplier1Entry = self.getWidget("multiplier1_label")
	#	self.update()


	def update( self ):
        
		aString = str( self.theFullPN()[ID] )
		aString += ':\n' + str( self.theFullPN()[PROPERTY] )        
		self.theIDEntry.set_text  ( aString )

		aValue = self.theSession.theSimulator.getProperty( createFullPNString( self.theFullPN() ) )
		

		value = aValue
		self.theActualValue = value
		self.theBarLength , self.theMultiplier , self.thePositiveFlag \
		                              = self.calculateBarLength( value )

		aIndicator = (value / (float)(10**(self.theMultiplier))) \
		                              * self.thePositiveFlag

		self['progressbar'].set_value(int(self.theBarLength))
		self['progressbar'].set_format_string(str(value))

		self.theMultiplier1Entry.set_text(str(int(self.theMultiplier-1)))
		self['multiplier_entry'].set_text(str(int(self.theMultiplier+2)))


	def updateByAuto( self, value ):

		self.theAutoChangeFlag = 1
		self.update()


	def updateByAddbutton( self , obj ):

		self['auto_button'].set_active( 0 )
		aNumberString =  self['multiplier_entry'].get_text()
		aNumber = string.atof( aNumberString )
		aNumber = aNumber + 1
		self.pull = aNumber

		self.theAutoChangeFlag = 0
		self.update()


	def updateBySubtractbutton( self,obj ):

		self['auto_button'].set_active( 0 )
		aNumberString =  self['multiplier_entry'].get_text()
		aNumber = string.atof( aNumberString )
		aNumber = aNumber - 1
		self.pull = aNumber
		self.theAutoChangeFlag = 0
		self.update()


	def updateByTextentry(self, obj):

		if self.theAutoChangeFlag :
			pass
		else :
			self['auto_button'].set_active( 0 )

		aNumberString = obj.get_text()

		aNumber = string.atof( aNumberString )
		self.pull = aNumber
		self.theAutoChangeFlag = 0
		self.update()


	def updateByAutoButton(self, autobutton):
		self.update()


	def calculateBarLength( self, value ):

		if value < 0 :
			value = - value
			aPositiveFlag = -1
		else :
			aPositiveFlag = 1

		if self['auto_button'].get_active() :
			if value == 0 :
				aMultiplier = 2
			else :
				aMultiplier = (int)(Numeric.log10(value))
			self.pull = aMultiplier+2
		else :
			aMultiplier = self.pull-2

		if value == 0:
			aBarLength = 0
		else :
			aBarLength = (Numeric.log10(value)+1-aMultiplier)*1000/3

		return  aBarLength, aMultiplier, aPositiveFlag
                

	def changeValue( self, value ):
		self.updateByAuto( value )

   
if __name__ == "__main__":

    class simulator:

        dic={('Substance','/CELL/CYTOPLASM','ATP','quantity') : (1950,),}
        
        def getProperty( self, fpn ):
            return simulator.dic[fpn]
        
        def setProperty( self, fpn, value ):
            simulator.dic[fpn] = value


    fpn = ('Substance','/CELL/CYTOPLASM','ATP','quantity')

    def mainQuit( obj, data ):
        gtk.mainquit()

    def mainLoop():
        # FIXME: should be a custom function
        gtk.mainloop()

    def main():
        aPluginManager = Plugin.PluginManager()
        aBargraphWindow = BargraphWindow( 'plugins', simulator(), [fpn,], aPluginManager )
        aBargraphWindow.addHandler( 'gtk_main_quit', mainQuit )
        
        mainLoop()


    # propertyValue1 = -750.0000

    main()
