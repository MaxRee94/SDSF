import sys
import time
import json

from widgets import SettingsWidget
from gui_utils import HContainer, VContainer, ListSettingWidget, KEYCODES
from PySide6 import QtCore, QtWidgets, QtGui

import config

module = sys.modules[__name__]
module.window = None

class MinimalSettingsWidget(SettingsWidget):
    """Slightly modified SettingsWidget that works without staging-/browsing widgets.
    """
    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self.settings_container.layout.removeRow(0)
        self.settings_container.path_widget = QtWidgets.QLineEdit()
        self.settings_container.path_widget.setPlaceholderText("Enter publish path...")
        self.settings_container.layout.insertRow(0, self.settings_container.path_widget)


class MinimalWindow(QtWidgets.QDialog):
    """Main window of minimal GUI."""

    def __init__(self, parent=None):
        super(MinimalWindow, self).__init__(parent)
        self.setWindowTitle("DispRes")

        # Enable minimize and maximize for app
        self.setWindowFlags(QtCore.Qt.Window)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)

        self.settings_widget = MinimalSettingsWidget()

        start_container = HContainer()
        start_btn = QtWidgets.QPushButton("Start")
        start_btn.setFixedWidth(120)
        start_container.layout.addStretch()
        start_container.layout.addWidget(start_btn)
        start_container.layout.addStretch()
        self.status_label = QtWidgets.QLabel("")

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.settings_widget)
        layout.addWidget(start_container)
        layout.addWidget(self.status_label)

        # Connections
        start_btn.clicked.connect(self._on_start_clicked)

        # Defaults
        self.resize(1000, 650)
        self.load_defaults()

    def refresh(self):
        self.echo("Refreshing..")

    def echo(self, message):
        print(message)

    def load_defaults(self):
        print("setting defaults..")

    def set_default_families(self):
        families_listwidget = (
            self.settings_widget.settings_container.families_listwidget
        )
        for fam in config.gui_defaults.get("families"):
            families_listwidget.listwidget.insertItem(0, fam)

    def set_default_family(self):
        self.settings_widget.settings_container.family_widget.setText(
            config.gui_defaults.get("family")
        )

    def collect_settings(self):
        """Collect and return publish settings entered by user."""
        path = self.settings_widget.settings_container.path_widget.text() or None
        family = self.settings_widget.settings_container.family_widget.text() or None
        families = self.settings_widget.settings_container.families_listwidget.get_items_texts() or None
        source = self.settings_widget.settings_container.source_line_edit.text() or None
        custom = self.settings_widget.settings_container.custom_instance_data_widget.toPlainText() or None

        settings = {
            "gridsize": [path],
            "family": family,
            "families": families,
            "source": source,
            "additional_instance_data": custom
        }

        return settings

    def _on_start_clicked(self):
        self.status_label.setText("Publishing...")

        print("Preparing publish..")
        settings = self.collect_settings()
        print("Publish Settings:", settings)
        try:
            #propeller.app.main(**settings)
            self.status_label.setText("Publish Succesful.")
        except RuntimeError:
            self.status_label.setText(
                "Publish failed. View terminal for debug messages"
            )


def show(debug=False):
    """Display Minimal GUI

    Args:
        debug (bool, optional): Run GUI in debug-mode, defaults to False
    """
    try:
        module.window.close()
        del module.window
    except (RuntimeError, AttributeError):
        pass

    
    window = MinimalWindow()
    window.show()

    window.refresh()

    module.window = window

    # Pull window to the front.
    module.window.raise_()
    module.window.activateWindow()

    return module.window
    

if __name__ == "__main__":
    app = QtWidgets.QApplication()
    window = show()
    app.exec()
    app.deleteLater()
