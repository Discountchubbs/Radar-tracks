/****************************************************************************
** Meta object code from reading C++ file 'map.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.12.8)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "map.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'map.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.12.8. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_MapWindow_t {
    QByteArrayData data[11];
    char stringdata0[107];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MapWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MapWindow_t qt_meta_stringdata_MapWindow = {
    {
QT_MOC_LITERAL(0, 0, 9), // "MapWindow"
QT_MOC_LITERAL(1, 10, 14), // "requestVisible"
QT_MOC_LITERAL(2, 25, 0), // ""
QT_MOC_LITERAL(3, 26, 7), // "visible"
QT_MOC_LITERAL(4, 34, 15), // "requestSettings"
QT_MOC_LITERAL(5, 50, 8), // "settings"
QT_MOC_LITERAL(6, 59, 11), // "updateState"
QT_MOC_LITERAL(7, 71, 7), // "UIState"
QT_MOC_LITERAL(8, 79, 1), // "s"
QT_MOC_LITERAL(9, 81, 17), // "offroadTransition"
QT_MOC_LITERAL(10, 99, 7) // "offroad"

    },
    "MapWindow\0requestVisible\0\0visible\0"
    "requestSettings\0settings\0updateState\0"
    "UIState\0s\0offroadTransition\0offroad"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MapWindow[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   34,    2, 0x06 /* Public */,
       4,    1,   37,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       6,    1,   40,    2, 0x08 /* Private */,
       9,    1,   43,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Bool,    3,
    QMetaType::Void, QMetaType::Bool,    5,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 7,    8,
    QMetaType::Void, QMetaType::Bool,   10,

       0        // eod
};

void MapWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<MapWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->requestVisible((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 1: _t->requestSettings((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 2: _t->updateState((*reinterpret_cast< const UIState(*)>(_a[1]))); break;
        case 3: _t->offroadTransition((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (MapWindow::*)(bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MapWindow::requestVisible)) {
                *result = 0;
                return;
            }
        }
        {
            using _t = void (MapWindow::*)(bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&MapWindow::requestSettings)) {
                *result = 1;
                return;
            }
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject MapWindow::staticMetaObject = { {
    &QOpenGLWidget::staticMetaObject,
    qt_meta_stringdata_MapWindow.data,
    qt_meta_data_MapWindow,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *MapWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MapWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MapWindow.stringdata0))
        return static_cast<void*>(this);
    return QOpenGLWidget::qt_metacast(_clname);
}

int MapWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QOpenGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 4)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 4)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void MapWindow::requestVisible(bool _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void MapWindow::requestSettings(bool _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
