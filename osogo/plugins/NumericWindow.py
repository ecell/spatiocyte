#!/usr/bin/env python


from PluginWindow import *
from ecssupport import *
import GTK

class NumericWindow( PluginWindow ):

    def __init__( self, dirname, sim, data, pluginmanager ):

        PluginWindow.__init__( self, dirname, sim, data, pluginmanager )

        self['toolbar5'].set_style( GTK.TOOLBAR_ICONS )
        self['toolbar5'].set_button_relief( GTK.RELIEF_HALF )

        self.addHandlers( { 'input_value'    :self.inputValue,
                            'increase_value' :self.increaseValue,
                            'decrease_value' :self.decreaseValue } )

        aPropertyListFullPN = \
                  convertFullIDToFullPN(self.theFullID(),'PropertyList')
        aPropertyList = \
                  list( self.theSimulator.getProperty( aPropertyListFullPN ) )

        aAttributeListFullPN = \
                  convertFullIDToFullPN(self.theFullID(), 'PropertyAttributes')
        aAttributeList = \
                  list(self.theSimulator.getProperty( aAttributeListFullPN ))

        num = 0
        for aProperty in aPropertyList:
            if (aProperty == 'Quantity' ):
                print aProperty,
                print "=",
                print aAttributeList[num]
            else :
                pass
            num += 1        

        self.initialize( self.theFullPN() )

    def initialize( self, fpn ):
        aString = str( self.theFullPN()[ID] )
        aString += ':\n' + str( self.theFullPN()[PROPERTY] )
        self["id_label"].set_text( aString )
        self.update()
        
    def update( self ):
        self["value_frame"].set_text( str( self.getValue() ) )

    def inputValue( self, obj ):
        aValue =  string.atof( obj.get_text() )
        self.setValue( aValue )

    def increaseValue( self, obj ):
        self.setValue( self.getValue() * 2.0 )

    def decreaseValue( self, obj ):
        self.setValue( self.getValue() * 0.5 )

    def getValue( self ):
        aValueList = self.theSimulator.getProperty( self.theFullPN() )
        return aValueList[0]

    def setValue( self, aValue ):
        aValueList = ( aValue, )
        self.theSimulator.setProperty( self.theFullPN(), aValueList )
        self.thePluginManager.updateAllPluginWindow()

### test code

if __name__ == "__main__":

    class simulator:

        dic={('Substance', '/CELL/CYTOPLASM', 'ATP','Quantity') : (1950,),}

        def getProperty( self, fpn ):
            return simulator.dic[fpn]

        def setProperty( self, fpn, value ):
            simulator.dic[fpn] = value


    fpn = ('Substance','/CELL/CYTOPLASM','ATP','Quantity')

    def mainQuit( obj, data ):
        print obj,data
        gtk.mainquit()
        
    def mainLoop():
        # FIXME: should be a custom function

        gtk.mainloop()

    def main():
        aNumericWindow = NumericWindow( 'plugins', simulator(), [fpn,] )
        aNumericWindow.addHandler( 'gtk_main_quit', mainQuit )
        aNumericWindow.update()

        mainLoop()

    main()









